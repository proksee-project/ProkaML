'''
Copyright:

University of Manitoba & National Microbiology Laboratory, Canada, 2020

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
import sys
from Bio import Entrez
from get_genomic_metadata import GetMetadata

if len(sys.argv) != 6:
    sys.exit('''
    Command usage: python metadata_extract_fileindex.py EMAIL NCBI_API_KEY INPUT_DIRECTORY OUTPUT_DIRECTORY FILE_NAME.
    Need to pass 5 arguments corresponding to your email, ncbi api key, input directory of files containing
    UIDs, output dirctory and a file name containing list of UIDs.
    ''')

else:
    email = sys.argv[1]
    api_key = sys.argv[2]
    input_dir = sys.argv[3]
    output_dir = sys.argv[4]
    filename = sys.argv[5]

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(os.path.join(input_dir, filename)) as f:
        idlist = f.read().splitlines()

    output_filename = filename.split('idlist.txt')[0] + 'metadata.txt'
    output_file = open(os.path.join(output_dir, output_filename), 'w')

    Entrez.email = email
    Entrez.api_key = api_key
    idlist_metadata = GetMetadata(idlist)
    idlist_metadata.print_genomic_metadata(output_file)
