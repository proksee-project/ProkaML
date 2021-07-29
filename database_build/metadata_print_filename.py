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
from get_genomic_metadata import AttributeMetadata

ID_LIST_FILE_EXTENSION = 'idlist.txt'
METADATA_FILE_EXTENSION = 'metadata.txt'

my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key input_file_path output_dir',
                                    description='Obtains assembly attributes from API queries')
my_parser.add_argument('email',
                        type=str,
                        help='user email address')
my_parser.add_argument('api_key',
                        type=str,
                        help='NCBI user API key')
my_parser.add_argument('input_file_path',
                        type=str,
                        help='path to file containing UIDs')
my_parser.add_argument('output_dir',
                        type=str,
                        help='output directory')
args = my_parser.parse_args()

input_file_path = args.input_file_path
output_dir = args.output_dir

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

with open(input_file_path) as f:
    idlist = f.read().splitlines()

output_filename = os.path.basename(input_file_path).split(ID_LIST_FILE_EXTENSION)[0] + METADATA_FILE_EXTENSION
output_file = open(os.path.join(output_dir, output_filename), 'w')

Entrez.email = args.email
Entrez.api_key = args.api_key
idlist_metadata = AttributeMetadata(idlist)
idlist_metadata.print_genomic_metadata(output_file)
