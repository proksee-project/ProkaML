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

import pandas as pd
from Bio import Entrez
from datetime import date
import argparse
import re

Entrez.max_tries = 10
INPUT_FILE_PREFIX = 'species_assemblycounts_'
FILE_EXTENSION = '.txt'
OUTPUT_FILE_SUFFIX = '_taxonomy'
SEPARATOR = '\t'


def get_taxonomy_kingdom(email, api_key, species):
    Entrez.email = email
    Entrez.api_key = api_key

    try:
        handle = Entrez.esearch(db="Taxonomy", term=species)
        record = Entrez.read(handle)

        efetch = Entrez.efetch(db="Taxonomy", id=record['IdList'], retmode="xml")
        taxon_summary = Entrez.read(efetch)

        taxonomic_dict = {}
        for i in range(len(taxon_summary[0]['LineageEx'])):
            taxonomic_dict[taxon_summary[0]['LineageEx'][i]['Rank']] = taxon_summary[0]['LineageEx'][i]['ScientificName']

        if 'superkingdom' in taxonomic_dict:
            kingdom = taxonomic_dict['superkingdom']

        else:
            kingdom = 'NA (Not available)'

    except Exception:
        kingdom = 'NA (Cannot be obtained)'

    return kingdom


def main():
    my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key',
                                        description='Gets full taxonomical lineage for species')
    my_parser.add_argument('email',
                           type=str,
                           help='user email address')
    my_parser.add_argument('api_key',
                           type=str,
                           help='NCBI user API key')                      
    args = my_parser.parse_args()

    email = args.email
    api_key = args.api_key

    month_year_stamp = date.today().strftime("%b_%Y")
    input_file = INPUT_FILE_PREFIX + month_year_stamp + FILE_EXTENSION
    dataframe = pd.read_csv(input_file, sep=SEPARATOR)
    species_list = dataframe['Species'].to_list()
    num_assemblies = dataframe['Num_assemblies'].to_list()
    output_file = open(INPUT_FILE_PREFIX + month_year_stamp + OUTPUT_FILE_SUFFIX + FILE_EXTENSION, 'w')
    output_file.write('Species\tNum_assemblies\tKingdom\n')
    for i in range(len(species_list)):
        species_corrected = re.sub(r'(\[|\]|\(.+?\))', '', species_list[i])
        kingdom = get_taxonomy_kingdom(email, api_key, species_corrected)
        print(str(i) + " species' kingdom retrieved")
        output_file.write('{}\t{}\t{}\n'.format(species_corrected ,num_assemblies[i], kingdom))


if __name__ == '__main__':
	main()
