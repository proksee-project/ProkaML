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

from pathlib import Path, PurePath
parent_directory = Path(__file__).resolve().parents[1]
ADDITIONAL_GENOMIC_ATTRIBUTES_PATH = PurePath.joinpath(parent_directory, 'add_genomic_attributes')

import sys
sys.path.append(ADDITIONAL_GENOMIC_ATTRIBUTES_PATH)

import os
import argparse
import pandas as pd
from entrez_metadata import EntrezMetadata
from gc_content import GCContentCalculate

SEPARATOR = '\t'
FILENAME_SPLIT_PATTERN = '_metadata.txt'
FILENAME_ID_INDEX = 0
FILENAME_EXTN = '_added_attributes.txt'
GC_CONTENT = 'GCcontent'
LOW_MEMORY = False
KEEP_DEFAULT_NA = False
WRITE_MODE = 'w'
KEEP_INDEX = False

ADDITIONAL_METADATA_DIR = 'additional_species_metadata'
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
species_name = input_file.split(FILENAME_SPLIT_PATTERN)[FILENAME_ID_INDEX]
output_file = input_file.split(FILENAME_SPLIT_PATTERN)[FILENAME_ID_INDEX] + FILENAME_EXTN
output_file_path = os.path.join(ADDITIONAL_METADATA_DIR, output_file)
species_log_file = open(species_name + '_log_gc_content.txt', WRITE_MODE)

dataframe = pd.read_csv(input_file_path, sep=SEPARATOR, keep_default_na=False, low_memory=LOW_MEMORY)
calculate_gc_content = GCContentCalculate(email, api_key, dataframe, species_log_file)
dataframe[GC_CONTENT], num_success, num_failure = calculate_gc_content.append_gc_content()

dataframe.to_csv(output_file_path, sep=SEPARATOR, mode=WRITE_MODE, index=KEEP_INDEX)

species_log_file.write('{}\t{}\n'.format(num_success, num_failure))
