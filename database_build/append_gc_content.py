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
import pandas as pd
import constants as const
from entrez_metadata import EntrezMetadataDownloader
from manage_assembly import AssemblyManager


def main():
    """
    Appends column of GC content of assembly
    """

    filename_split_pattern = const.Assembly.METADATA_SUFFIX + const.FileFormat.TEXT
    FILENAME_ID_INDEX = 0
    LOW_MEMORY = False
    KEEP_DEFAULT_NA = False
    KEEP_INDEX = False

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
                            help='path to file containing species metadata')

    args = my_parser.parse_args()

    email = args.email
    api_key = args.api_key
    input_file_path = args.input_file_path

    input_file = os.path.basename(input_file_path)
    species_name = input_file.split(filename_split_pattern)[FILENAME_ID_INDEX]
    output_file = input_file.split(filename_split_pattern)[FILENAME_ID_INDEX] + const.Assembly.GC_SUFFIX + const.FileFormat.TEXT
    output_file_path = os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.GC_METADATA_DIR, output_file)
    species_log_file = open(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.SPECIES_GC_LOG_DIR,
        (species_name + const.LogFiles.SUB_LOG_GC + const.FileFormat.TEXT)), const.FileFormat.WRITE_MODE)

    dataframe = pd.read_csv(input_file_path, sep=const.FileFormat.SEPARATOR, keep_default_na=KEEP_DEFAULT_NA, low_memory=LOW_MEMORY)
    dataframe[const.Metadata.METADATA_COLUMN_HEADERS[const.Metadata.METADATA_INDEX_GC_CONTENT]], num_success, fail_server, \
        fail_mismatch, fail_zlib = AssemblyManager(email, api_key, dataframe, species_log_file).append_gc_content()
    dataframe.to_csv(output_file_path, sep=const.FileFormat.SEPARATOR , mode=const.FileFormat.WRITE_MODE, index=KEEP_INDEX)
    species_log_file.write('{}\t{}\t{}\t{}\n'.format(num_success, fail_server, fail_mismatch, fail_zlib))


if __name__ == '__main__':
    main()
