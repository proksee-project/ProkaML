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
import argparse
from Bio import Entrez
from entrez_metadata import EntrezMetadata
import constants as const
import re


def main():
    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key input_file_path',
                                        description='Obtains assembly attributes from API queries')
    my_parser.add_argument('email',
                            type=str,
                            help='user email address')
    my_parser.add_argument('api_key',
                            type=str,
                            help='NCBI user API key')
    my_parser.add_argument('input_file_path',
                            type=str,
                            help='path to file containing assembly UIDs')

    args = my_parser.parse_args()

    input_file_path = args.input_file_path
    file_count_pattern = re.search(r'.*chunk(\d+)\.txt', os.path.basename(input_file_path))
    file_number = file_count_pattern.group(1)
    log_file =  open(os.path.join(const.FileDirectories.DATABASE_PATH, 'log_entrez_metadata_chunk' + file_number + \
        const.FileFormat.FILE_EXTENSION), const.FileFormat.WRITE_MODE)

    with open(input_file_path) as f:
        id_list = f.read().splitlines()

    output_filename = os.path.basename(input_file_path).split('.')[0] + const.Assembly.METADATA_SUFFIX + const.FileFormat.FILE_EXTENSION
    output_file = open(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ENTREZ_METADATA_DIR, output_filename), \
        const.FileFormat.WRITE_MODE)

    Entrez.email = args.email
    Entrez.api_key = args.api_key
    idlist_metadata = EntrezMetadata(id_list, output_file)
    num_success_uids, irretrievable_uids = idlist_metadata.obtain_log_metadata()
    log_file.write('{}\t{}\t{}\n'.format(file_number, num_success_uids, irretrievable_uids))
    log_file.close()


if __name__ == '__main__':
    main()
