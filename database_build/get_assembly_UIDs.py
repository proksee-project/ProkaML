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
from datetime import date
import argparse
import constants as const
import re
import os


def count_overall_assemblydb_number(email, api_key, log):

    Entrez.email = email
    Entrez.api_key = api_key

    for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
        try:
            handle = Entrez.esearch(db=const.ASSEMBLY_DATABASE, term=const.CONTIG_TERM)

        except Exception:
            overall_num_assemblies = 0

        else:
            record = Entrez.read(handle)
            overall_num_assemblies = int(record[const.RECORD_COUNT])
            break

    if overall_num_assemblies == 0:
        log_message = 'Either internet failure OR invalid/incorrect ordering ' + \
            'of input parameters (email, API key) OR some unexpected failure..Exiting...\n'
        log.write(log_message)

    return overall_num_assemblies


def write_assembly_UIDs(overall_num_assemblies):

    if overall_num_assemblies > 0:
        for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
            try:
                handle1 = Entrez.esearch(db=const.ASSEMBLY_DATABASE, 
                                         term=const.CONTIG_TERM, 
                                         retmax=overall_num_assemblies)

            except Exception:
                num_batches = 0

            else:
                record1 = Entrez.read(handle1)
                num_batches, last_batch_size = divmod(overall_num_assemblies, const.ESUMMARY_BATCH_LIMIT)

                for batch_instance in range(num_batches + 1):
                    UID_output_file = open(os.path.join(const.UID_OUTPUT_DIR, const.UID_PREFIX + \
                        str(int(batch_instance + 1)) + const.FILE_EXTENSION), const.WRITE_MODE)
                    batch_start = const.ESUMMARY_BATCH_LIMIT * batch_instance

                    if batch_instance < num_batches:
                        batch_end = batch_start + const.ESUMMARY_BATCH_LIMIT

                    elif batch_instance == num_batches:
                        batch_end = batch_start + last_batch_size

                    for i in range(batch_start, batch_end):
                        UID_output_file.write(record1[const.RECORD_IDLIST][i] + const.NEW_LINE)

                break

    else:
        num_batches = 0

    if num_batches == 0:
        return num_batches
    else:
        return int(num_batches + 1)


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

    log = open(const.LOG_FILE + const.FILE_EXTENSION, mode=const.WRITE_MODE)
    log.write('#########################################################\n')
    log.write('Getting NCBI assembly UIDs in batches of 10,000\n')
    log.write('#########################################################\n')

    num_assemblies = count_overall_assemblydb_number(email, api_key, log)
    if num_assemblies > 0:
        if not os.path.exists(const.UID_OUTPUT_DIR):
            os.mkdir(const.UID_OUTPUT_DIR)
        if not os.path.exists(const.ENTREZ_METADATA_DIR):
            os.mkdir(const.ENTREZ_METADATA_DIR)
    num_batches = write_assembly_UIDs(num_assemblies)
    log.write('{} NCBI assembly UIDs written to {} file chunks - format {}/AssemblyUID_chunk*.txt\n'.format(num_assemblies, num_batches, const.UID_OUTPUT_DIR))
    log.close()


if __name__ == '__main__':
    main()
