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

from Bio import Entrez
from collections import defaultdict
import sys
from datetime import date
import argparse

DATABASE = 'assembly'
TERM = 'contig[Assembly Level]'
COUNT = 'Count'
BATCH_LIMIT = 10000
ID_LIST = 'IdList'
ID_LIST_JOIN_CHAR = ','
VALIDATE = False
DOCUMENT_SUMMARY_SET = 'DocumentSummarySet'
DOCUMENT_SUMMARY = 'DocumentSummary'
SPECIES = 'SpeciesName'
OUTPUT_FILE_PREFIX = 'species_assemblycounts_'
OUTPUT_FILE_EXTENSION = '.txt'
ASSEMBLY_LOWERBOUND = 10
ERROR_LOG_FILE = open('error_log_species_assemblies.txt', 'w')
ERROR_MESSAGE_SERVER = 'NCBI server connection error'


def count_assem_records(email, api_key):

    Entrez.email = email
    Entrez.api_key = api_key

    try:
        handle = Entrez.esearch(db=DATABASE, term=TERM)
        record = Entrez.read(handle)
        num_assemblies = record[COUNT]

    except Exception:
        raise Exception('Either internet failure OR invalid/incorrect ordering ' +
            'of input parameters (email, API key) OR some unexpected failure..Exiting...')

    return int(num_assemblies)


def retrieve_assem_records(num_assemblies):

    species_assembly_dict = defaultdict(int)

    handle1 = Entrez.esearch(db=DATABASE, term=TERM, retmax=num_assemblies)
    record1 = Entrez.read(handle1)
    num_batches, last_batch_size = divmod(num_assemblies, BATCH_LIMIT)

    for batch_instance in range(0, num_batches + 1):
        idlist_batch_instance = []
        batch_start = BATCH_LIMIT * batch_instance

        if batch_instance < num_batches:
            batch_end = batch_start + BATCH_LIMIT

        elif batch_instance == num_batches:
            batch_end = batch_start + last_batch_size

        for i in range(batch_start, batch_end):
            idlist_batch_instance.append(record1[ID_LIST][i])

        append_species_dict(idlist_batch_instance, species_assembly_dict)

        string1 = str(batch_end) + ' document summaries processed. '
        string2 = 'Species dictionary has ' + str(len(species_assembly_dict)) + ' records'
        print(string1 + string2)

    return species_assembly_dict


def append_species_dict(idlist_batch_instance, species_assembly_dict):

    esum = Entrez.esummary(db=DATABASE, id=ID_LIST_JOIN_CHAR.join(idlist_batch_instance))
    docsum = Entrez.read(esum, validate=VALIDATE)

    for j in range(0, len(idlist_batch_instance)):

        try:
            species = docsum[DOCUMENT_SUMMARY_SET][DOCUMENT_SUMMARY][j][SPECIES]
            species_assembly_dict[species] += 1

        except IndexError:
            # Species can't be retrieved due to NCBI server connection issue. Skipping and moving ahead
            error = genbank_id + ' ' + ERROR_MESSAGE_SERVER
            ERROR_LOG_FILE.write(error)

            pass


def species_dicn_write(species_assembly_dict, output_file):

    if species_assembly_dict:
        for species in sorted(species_assembly_dict, key=species_assembly_dict.get, reverse=True):
            if species_assembly_dict[species] >= 10:
                output_file.write('{}\t{}\n'.format(species, species_assembly_dict[species]))

        success = 'Species dictionary written'

    else:
        success = 'No species dictionary from previous steps'

    return success


def main():

    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key',
                                        description='Runs Entrez API queries on NCBI assembly database')
    my_parser.add_argument('email',
                           type=str,
                           help='user email address')
    my_parser.add_argument('api_key',
                           type=str,
                           help='NCBI user API key')                      
    args = my_parser.parse_args()

    email = args.email
    api_key = args.api_key

    try:
        num_assemblies = count_assem_records(email, api_key)
    except Exception as e:
        sys.exit(e)

    species_assembly_dict = retrieve_assem_records(num_assemblies)

    month_year_stamp = date.today().strftime("%b_%Y")
    output_file = open(OUTPUT_FILE_PREFIX + month_year_stamp + OUTPUT_FILE_EXTENSION, 'w')
    output_file.write('Species\tNum_assemblies\n')

    print(species_dicn_write(species_assembly_dict, output_file))


if __name__ == '__main__':
    main()
