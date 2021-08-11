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

CHUNK_DIVISOR = 10000
SEPARATOR = '\t'
KINGDOM = 'Kingdom'
SPECIES = 'Species'
NUM_ASSEMBLIES = 'Num_assemblies'
KINGDOM_FILTER = 'Bacteria'

DATABASE = 'assembly'
TERM_EXTENSION = '[Organism] AND contig[Assembly Level]'
ID_LIST = 'IdList'
ID_LIST_LIMIT = 10000
OUTPUT_DIR = 'id_list'


def select_bacterial_species(input_file):

    # Subsets dataframe for assemblies belonging to bacterial kingdom
    dataframe = pd.read_csv(input_file, sep=SEPARATOR)
    bacterial_dataframe = dataframe[dataframe[KINGDOM] == KINGDOM_FILTER]

    return bacterial_dataframe


def get_id_list(species, num_assemblies):

    # Returns list of assembly UIDs for a given species
    try:
        handle = Entrez.esearch(db=DATABASE,
                                term=species + TERM_EXTENSION,
                                retmax=num_assemblies)
        record = Entrez.read(handle)

    except Exception:
        raise Exception('Either connection failure OR some unexpected failure..Exiting...')

    return record[ID_LIST]


def get_num_chunks(id_list):

    # Calculates number of file chunks for a species' list of UIDs
    if len(id_list) % ID_LIST_LIMIT == 0:
        num_chunks = int(len(id_list)/ID_LIST_LIMIT)
    else:
        num_chunks = int(len(id_list)/ID_LIST_LIMIT) + 1

    return num_chunks


def get_species_idlist_filelist(species, num_chunks):

    # Generates species-chunk specific output files for writing UIDs
    species_outfile_list = []
    species_name_list = re.split(r'\s', species)
    for i in range(1, num_chunks + 1):
        outfile = "_".join(species_name_list) + '_chunk' + str(i) + '_idlist.txt'
        species_outfile_list.append(outfile)

    return species_outfile_list


def write_idlist_filechunks(index, species, id_list, num_chunks, file_list):

    # Writes UIDs to species-chunk specific output files
    count = index + 1
    if num_chunks == 0:
        print(str(count) + ' species processing: ' + str(species) +
                ' returns empty ID list with esearch. skipping..')

    elif num_chunks > 0:
        print(str(count) + ' species processing: ' + str(species))
        for k in range(0, num_chunks):
            output_file = open(os.path.join(OUTPUT_DIR, file_list[k]), 'w')
            idlist_lowerlimit = k * ID_LIST_LIMIT
            idlist_upperlimit = idlist_lowerlimit + ID_LIST_LIMIT
            for m in range(idlist_lowerlimit, min(idlist_upperlimit, len(id_list))):
                output_file.write('{}\n'.format(id_list[m]))


def main():
    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key input_file',
                                        description='Retrieves assembly UIDs from API queries')
    my_parser.add_argument('email',
                           type=str,
                           help='user email address')
    my_parser.add_argument('api_key',
                           type=str,
                           help='NCBI user API key')
    my_parser.add_argument('input_file',
                           type=str,
                           help='Input file containing species names, assembly counts and ' +
                                'taxonomical kingdom')
    args = my_parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = args.api_key
    input_file = args.input_file
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    bacterial_dataframe = select_bacterial_species(input_file)
    species_list = bacterial_dataframe[SPECIES].to_list()
    num_assemblies_list = bacterial_dataframe[NUM_ASSEMBLIES].to_list()

    for i in range(0, len(species_list)):
        species_id_list = get_id_list(species_list[i], num_assemblies_list[i])
        species_num_chunks = get_num_chunks(species_id_list)
        species_outfile_list = get_species_idlist_filelist(species_list[i], species_num_chunks)

        write_idlist_filechunks(i, species_list[i], species_id_list, species_num_chunks, species_outfile_list)


if __name__ == "__main__":
    main()
