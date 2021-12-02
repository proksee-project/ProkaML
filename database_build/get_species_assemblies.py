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

Entrez.max_tries = 1
Entrez.sleep_between_tries = 1
ASSEMBLY_DATABASE = 'Assembly'
TERM = 'contig[Assembly Level]'
COUNT = 'Count'
BATCH_LIMIT = 10000
ID_LIST = 'IdList'
ID_LIST_JOIN_CHAR = ','
ESUMMARY_VALIDATE = False
DOCUMENT_SUMMARY_SET = 'DocumentSummarySet'
DOCUMENT_SUMMARY = 'DocumentSummary'
SPECIES = 'SpeciesName'
OUTPUT_FILE_PREFIX = 'species_assemblycounts_'
OUTPUT_FILE_EXTENSION = '.txt'
ASSEMBLY_COUNT_LOWERBOUND = 10
LOG_FILE = open('LOG.txt', 'w')


def count_overall_assemblydb_number(email, api_key):

    Entrez.email = email
    Entrez.api_key = api_key

    for attempts in range(1, 4):
        try:
            handle = Entrez.esearch(db=ASSEMBLY_DATABASE, term=TERM)

        except Exception:
            overall_num_assemblies = 0

        else:
            record = Entrez.read(handle)
            overall_num_assemblies = int(record[COUNT])
            if overall_num_assemblies > 0:
                break

    if overall_num_assemblies == 0:
        log_message = 'Either internet failure OR invalid/incorrect ordering ' + \
            'of input parameters (email, API key) OR some unexpected failure..Exiting...\n'
        LOG_FILE.write(log_message)

    return overall_num_assemblies


def retrieve_species_assembly_counts(overall_num_assemblies):

    species_assembly_dict = defaultdict(int)

    if overall_num_assemblies > 0:
        for attempts in range(1, 4):
            try:
                handle1 = Entrez.esearch(db=ASSEMBLY_DATABASE, term=TERM, retmax=overall_num_assemblies)
        
            except Exception:
                pass

            else:
                record1 = Entrez.read(handle1)
                num_batches, last_batch_size = divmod(overall_num_assemblies, BATCH_LIMIT)

                for batch_instance in range(0, num_batches + 1):
                    idlist_batch_instance = []
                    batch_start = BATCH_LIMIT * batch_instance

                    if batch_instance < num_batches:
                        batch_end = batch_start + BATCH_LIMIT

                    elif batch_instance == num_batches:
                        batch_end = batch_start + last_batch_size

                    for i in range(batch_start, batch_end):
                        idlist_batch_instance.append(record1[ID_LIST][i])

                    append_species_dict(idlist_batch_instance, species_assembly_dict, batch_start, batch_end)

                if species_assembly_dict:
                    break

    return species_assembly_dict


def append_species_dict(idlist_batch_instance, species_assembly_dict, batch_start_index, batch_end_index):

    esummary = {}
    for attempts in range(1, 4):
        try:
            esum = Entrez.esummary(db=ASSEMBLY_DATABASE, id=ID_LIST_JOIN_CHAR.join(idlist_batch_instance))

        except Exception:
            pass

        else:
            esummary[attempts] = 'success'
            docsum = Entrez.read(esum, validate=ESUMMARY_VALIDATE)
            for j in range(0, len(idlist_batch_instance)):
                try:
                    species = docsum[DOCUMENT_SUMMARY_SET][DOCUMENT_SUMMARY][j][SPECIES]

                except IndexError:
                    log_message = 'Species for assembly UID ' + str(idlist_batch_instance[j]) + \
                        ' is absent\n'
                    LOG_FILE.write(log_message)

                else:
                    species_assembly_dict[species] += 1

            break

    if esummary:
        log_message = 'Assembly UIDs indices ' + str(batch_start_index) + ' to ' + \
            str(batch_end_index) + "'s document summaries processed. Species " + \
            'dictionary has ' + str(len(species_assembly_dict)) + ' unique species\n'
        print(log_message, end='')

    else:
        log_message = 'esummary for assembly UIDs ' + str(batch_start_index) + ' to ' + \
            str(batch_end_index) + "'s cannot be retrieved due to NCBI server connection error\n"
        LOG_FILE.write(log_message)
        print(log_message, end='')


def write_species_dict(species_assembly_dict):

    if species_assembly_dict:
        num_species_to_file = 0
        month_year_stamp = date.today().strftime("%b_%Y")
        output_file_name = OUTPUT_FILE_PREFIX + month_year_stamp + OUTPUT_FILE_EXTENSION
        with open(output_file_name, 'w') as output_file:
            output_file.write('Species\tNum_assemblies\n')
            for species in sorted(species_assembly_dict, key=species_assembly_dict.get, reverse=True):
                if species_assembly_dict[species] >= ASSEMBLY_COUNT_LOWERBOUND:
                    output_file.write('{}\t{}\n'.format(species, species_assembly_dict[species]))
                    num_species_to_file += 1

        output_file.close()
        log_message = str(num_species_to_file) + ' species in NCBI assembly database with assembly ' + \
            'counts written to ' + output_file_name + '\n'

    else:
        log_message = "Species' assembly counts cannot be determined\n"

    return log_message


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

    LOG_FILE.write('\n#########################################################\n')
    LOG_FILE.write("Getting species' assembly counts in NCBI\n")
    LOG_FILE.write('#########################################################\n')
    
    num_assemblies = count_overall_assemblydb_number(email, api_key)
    species_assembly_dict = retrieve_species_assembly_counts(num_assemblies)
    final_message = write_species_dict(species_assembly_dict)
    LOG_FILE.write(final_message)
    LOG_FILE.close()


if __name__ == '__main__':
    main()
