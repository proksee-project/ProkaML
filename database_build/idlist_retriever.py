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

import os
from Bio import Entrez
import re
import argparse

CHUNK_DIVISOR = 10000


def species_select(input_filename):

    # Groups species counts into four categories depending on assembly counts for respective species
    MAJOR_SPECIES_THRESHOLD = 1000
    LARGE_SPECIES_THRESHOLD = 100
    INTERMEDIATE_SPECIES_THRESHOLD = 10

    species_list_major = []
    species_list_large = []
    species_list_interm = []
    species_list_minor = []

    infile = open(input_filename, 'r')
    for line in infile:
        species_line = line.rstrip().split('\t')
        species_corrected = re.sub(r'[\[\]]', '', species_line[0])
        species_tuple = species_corrected, species_line[1]

        if int(species_line[1]) >= MAJOR_SPECIES_THRESHOLD:
            species_list_major.append(species_tuple)
        elif int(species_line[1]) < MAJOR_SPECIES_THRESHOLD and int(species_line[1]) >= LARGE_SPECIES_THRESHOLD:
            species_list_large.append(species_tuple)
        elif int(species_line[1]) < LARGE_SPECIES_THRESHOLD and int(species_line[1]) >= INTERMEDIATE_SPECIES_THRESHOLD:
            species_list_interm.append(species_tuple)
        elif int(species_line[1]) < INTERMEDIATE_SPECIES_THRESHOLD:
            species_list_minor.append(species_tuple)

    species_entire_tuple = species_list_major, species_list_large, species_list_interm, species_list_minor

    return species_entire_tuple


def id_list(species_tuple=None):

    # Returns list of assembly UIDs for a given species
    DATABASE = 'assembly'
    TERM_EXTENSION = '[Organism] AND contig[Assembly Level]'
    ID_LIST = 'IdList'

    try:
        handle = Entrez.esearch(db=DATABASE,
                                term=species_tuple[0] + TERM_EXTENSION,
                                retmax=species_tuple[1])
        record = Entrez.read(handle)

    except Exception:
        raise Exception('Either connection failure OR some unexpected failure..Exiting...')

    return record[ID_LIST]


def chunk_counter(idlist=None):

    # Calculates number of file chunks for a species' list of UIDs
    CHUNK_DIVISOR = 10000

    if len(idlist) % CHUNK_DIVISOR == 0:
        chunks = len(idlist)/CHUNK_DIVISOR
    else:
        chunks = int(len(idlist)/CHUNK_DIVISOR) + 1

    return chunks


def species_outfile(species_tuple=None, chunks=None):

    # Generates chunk-wise output files for writing UIDs'''
    species_outfile_list = []
    species_name_list = re.split(r'\s|\/', species_tuple[0])
    for i in range(1, chunks + 1):
        outfile = "_".join(species_name_list) + '_chunk' + str(i) + '_idlist.txt'
        species_outfile_list.append(outfile)

    return species_outfile_list


def species_iterate(species_list=None, out_dir=None):

    '''Writes UIDs to species specific files'''
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for count, sp_tup in enumerate(species_list):
        idlist = id_list(sp_tup)
        chunks = chunk_counter(idlist)
        if chunks > 0:
            species_outfile_list = species_outfile(sp_tup, chunks)
            print(str(count) + ' species processing: ' + str(sp_tup[0]))
            for k in range(0, chunks):
                idlist_lowerlim = k*CHUNK_DIVISOR
                idlist_upperlim = (k+1)*CHUNK_DIVISOR
                idlist_batch = []
                for m in range(idlist_lowerlim, min(idlist_upperlim, len(idlist))):
                    idlist_batch.append(idlist[m])

                output_handle = os.path.join(out_dir, species_outfile_list[k])
                output_file = open(output_handle, 'w')
                idlist_string = '\n'.join(idlist_batch)
                output_file.write(idlist_string + '\n')

        elif chunks == 0:
            print(str(count) + ' species processing: ' + str(sp_tup[0]) +
                  ' returns empty ID list with esearch. skipping..')


def main():
    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key',
                                        description='Retrieves assembly UIDs from API queries')
    my_parser.add_argument('email',
                           type=str,
                           help='user email address')
    my_parser.add_argument('api_key',
                           type=str,
                           help='NCBI user API key')
    my_parser.add_argument('input_file',
                           type=str,
                           help='Input file containing two column tab separated entries for ' +
                                'species and assembly counts')
    args = my_parser.parse_args()
    
    Entrez.email = args.email
    Entrez.api_key = args.api_key
    input_file = args.input_file

    species_all_tuple = species_select(input_file)
    species_iterate(species_all_tuple[0], 'id_list_major')
    species_iterate(species_all_tuple[1], 'id_list_large')
    species_iterate(species_all_tuple[2], 'id_list_interm')
    species_iterate(species_all_tuple[3], 'id_list_minor')


if __name__ == "__main__":
    main()
