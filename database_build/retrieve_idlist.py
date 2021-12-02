'''
Copyright:

University of Manitoba & National Microbiology Laboratory, Canada, 2021

Written by: Arnab Saha Mandal

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
'''

import pandas as pd
import os
from Bio import Entrez
import re
import argparse
from datetime import date

INPUT_FILE_PREFIX = 'species_assemblycounts_'
INPUT_FILE_SUFFIX = '_taxonomy.txt'
OUTPUT_TABLE_SUFFIX = '_UIDs_numbers.txt'
SEPARATOR = '\t'
KINGDOM = 'Kingdom'
SPECIES = 'Species'
NUM_ASSEMBLIES = 'Num_assemblies'
NUM_UIDs = 'Num_UIDs'
PROKARYOTES = ['Bacteria', 'Archaea']

ASSEMBLY_DATABASE = 'Assembly'
TERM_EXTENSION = '[Organism] AND contig[Assembly Level]'
ID_LIST = 'IdList'
ID_LIST_LIMIT = 10000
OUTPUT_DIR = 'entrez_id_list'
ENTREZ_METADATA_DIR = 'entrez_species_metadata'
ADDITIONAL_METADATA_DIR = 'additional_species_metadata'
LOG_FILE = open('LOG.txt', 'a')


def select_prokaryote_species(input_file):

    # Subsets dataframe for assemblies belonging to prokaryotes
    dataframe = pd.read_csv(input_file, sep=SEPARATOR)
    prokaryote_dataframe = dataframe[dataframe[KINGDOM].isin(PROKARYOTES)]

    return prokaryote_dataframe


def get_id_list(species, num_assemblies):

    # Returns list of assembly UIDs for a given species
    for attempts in range(1, 4):
        try:
            handle = Entrez.esearch(db=ASSEMBLY_DATABASE,
                                    term=species + TERM_EXTENSION,
                                    retmax=num_assemblies)

        except Exception:
            record[ID_LIST] = []

        else:
            record = Entrez.read(handle)
            break

    return record[ID_LIST]


def get_species_num_chunks(species_id_list):

    # Calculates number of file chunks for a species' list of UIDs
    if len(species_id_list) % ID_LIST_LIMIT == 0:
        num_chunks = int(len(species_id_list)/ID_LIST_LIMIT)
    else:
        num_chunks = int(len(species_id_list)/ID_LIST_LIMIT) + 1

    return num_chunks


def get_species_idlist_filelist(species, num_chunks):

    # Generates species specific id list files with limit of 10000 id per file
    species_id_outfile_list = []
    species_split = re.split(r'\s', species)
    for i in range(1, num_chunks + 1):
        id_outfile = "_".join(species_split) + '_chunk' + str(i) + '_idlist.txt'
        species_id_outfile_list.append(id_outfile)

    return species_id_outfile_list


def write_idlist_filechunks(index, species, id_list, num_chunks, file_list):

    # Writes UIDs to species-chunk specific output files
    count = index + 1
    if num_chunks == 0:
        log_message = str(count) + ' species processing: ' + str(species) + \
            ' returns empty ID list\n'
        print(log_message, end='')

    elif num_chunks > 0:
        log_message = str(count) + ' species processing: ' + str(species) + \
            ' writing UIDs to file\n'
        print(log_message, end='')
        for k in range(0, num_chunks):
            output_file = open(os.path.join(OUTPUT_DIR, file_list[k]), 'w')
            idlist_lowerlimit = k * ID_LIST_LIMIT
            idlist_upperlimit = idlist_lowerlimit + ID_LIST_LIMIT
            for m in range(idlist_lowerlimit, min(idlist_upperlimit, len(id_list))):
                output_file.write('{}\n'.format(id_list[m]))


def main():
    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key',
                                        description='Retrieves assembly UIDs from API queries')
    my_parser.add_argument('email',
                           type=str,
                           help='user email address')
    my_parser.add_argument('api_key',
                           type=str,
                           help='NCBI user API key')
    args = my_parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = args.api_key

    month_year_stamp = date.today().strftime("%b_%Y")
    input_file = INPUT_FILE_PREFIX + month_year_stamp + INPUT_FILE_SUFFIX

    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    LOG_FILE.write('\n#########################################################\n')
    LOG_FILE.write("Retrieving Species' specific assembly UIDs\n")
    LOG_FILE.write('#########################################################\n')

    prokaryote_dataframe = select_prokaryote_species(input_file)
    log_message = str(prokaryote_dataframe.shape[0]) + ' species are prokaryotes (Bacteria/' + \
        'Archaea)\n' + 'Total assembly counts = ' + str(prokaryote_dataframe[NUM_ASSEMBLIES].sum()) + '\n'
    LOG_FILE.write(log_message)

    species_list = prokaryote_dataframe[SPECIES].to_list()
    num_assemblies_list = prokaryote_dataframe[NUM_ASSEMBLIES].to_list()
    num_UIDs_list = []

    total_difference = 0
    total_UIDs = 0
    for i in range(0, len(species_list)):
        species_id_list = get_id_list(species_list[i], num_assemblies_list[i])
        num_UIDs_list.append(len(species_id_list))

        difference = num_assemblies_list[i] - len(species_id_list)
        total_difference += difference
        total_UIDs += len(species_id_list)

        if difference > 0:
            log_message = 'Species ' + species_list[i] + ': ' + str(difference) + \
                ' UIDs could not be retrieved\n'
            LOG_FILE.write(log_message)

        species_num_chunks = get_species_num_chunks(species_id_list)
        species_id_outfile_list = get_species_idlist_filelist(species_list[i], species_num_chunks)
        write_idlist_filechunks(i, species_list[i], species_id_list, species_num_chunks, \
            species_id_outfile_list)

    log_message = 'Total UIDs retrieved = ' + str(total_UIDs) + '\nTotal UIDs not retrieved = ' \
        + str(total_difference) + '\n'
    LOG_FILE.write(log_message)

    if not os.path.exists(ENTREZ_METADATA_DIR):
        os.mkdir(ENTREZ_METADATA_DIR)
    if not os.path.exists(ADDITIONAL_METADATA_DIR):
        os.mkdir(ADDITIONAL_METADATA_DIR)


if __name__ == "__main__":
    main()
