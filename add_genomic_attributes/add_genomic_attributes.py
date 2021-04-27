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

import sys
from gc_content import GCContentCalculate
import argparse

"""
This program takes in a file (e.g. well_represented_species_metadata.txt) containing metadata of assembly
attributes, computes genomic attributes of interest for every assembly and appends those attributes
as additional columns within the file.
Currently, the only calculated attribute (in addition to NCBI derived attributes) is overall GC content,
which is the fraction of G or C bases of all nucleotide bases of an assembly.
"""

my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] EMAIL API_KEY FILE_NAME OUTPUT_DIRECTORY',
                                    description='Calculates genomic assembly attributes of interest')
my_parser.add_argument('email',
                        type=str,
                        help='user email address')
my_parser.add_argument('api_key',
                        type=str,
                        help='NCBI user API key')
my_parser.add_argument('filename',
                        type=str,
                        help='input file containing of assembly metadata')
my_parser.add_argument('output_dir',
                        type=str,
                        help='output directory')												
args = my_parser.parse_args()

email = args.email
api_key = args.api_key
filename = args.filename
output_dir = args.output_dir

calculate_gc_content = GCContentCalculate(email, api_key, filename, output_dir)
calculate_gc_content.append_gc_content()
