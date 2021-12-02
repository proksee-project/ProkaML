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

Entrez.max_tries = 1
Entrez.sleep_between_tries = 1
ID_LIST = 'IdList'
TAXONOMY_DATABASE = 'Taxonomy'
RETMODE = 'xml'
INPUT_FILE_PREFIX = 'species_assemblycounts_'
FILE_EXTENSION = '.txt'
OUTPUT_FILE_SUFFIX = '_taxonomy'
SEPARATOR = '\t'
LOG_FILE = open('LOG.txt', 'a')
LINEAGE = 'LineageEx'
RANK = 'Rank'
SCIENTIFIC_NAME = 'ScientificName'
SUPERKINGDOM = 'superkingdom'
PHYLUM = 'phylum'
CLASS = 'class'
ORDER = 'order'
FAMILY = 'family'
GENUS = 'genus'


def get_taxonomy_record(email, api_key, species):
    Entrez.email = email
    Entrez.api_key = api_key

    for attempts in range(1, 4):
        try:
            handle = Entrez.esearch(db=TAXONOMY_DATABASE, term=species)

        except Exception:
            record = {}

        else:
            record = Entrez.read(handle)
            if record:
                break

    return record


def get_full_taxonomy(record):

    taxonomy_dict = {}
    if record:
        for attempts in range(1, 4):        
            try:
                efetch = Entrez.efetch(db=TAXONOMY_DATABASE, id=record[ID_LIST], retmode=RETMODE)

            except Exception:
                pass
        
            else:
                taxon_summary = Entrez.read(efetch)
                for i in range(len(taxon_summary[0][LINEAGE])):
                    taxonomy_dict[taxon_summary[0][LINEAGE][i][RANK]] = taxon_summary[0][LINEAGE][i][SCIENTIFIC_NAME]

                break

    return taxonomy_dict


def get_taxonomy_kingdom(taxonomy_dict):

    if SUPERKINGDOM in taxonomy_dict:
        kingdom = taxonomy_dict[SUPERKINGDOM]

    else:
        kingdom = 'NA'

    return kingdom


def get_taxonomy_phylum(taxonomy_dict):

    if PHYLUM in taxonomy_dict:
        phylum = taxonomy_dict[PHYLUM]

    else:
        phylum = 'NA'

    return phylum


def get_taxonomy_class(taxonomy_dict):

    if CLASS in taxonomy_dict:
        taxonomy_class = taxonomy_dict[CLASS]

    else:
        taxonomy_class = 'NA'

    return taxonomy_class


def get_taxonomy_order(taxonomy_dict):

    if ORDER in taxonomy_dict:
        order = taxonomy_dict[ORDER]

    else:
        order = 'NA'

    return order


def get_taxonomy_family(taxonomy_dict):

    if FAMILY in taxonomy_dict:
        family = taxonomy_dict[FAMILY]

    else:
        family = 'NA'

    return family


def get_taxonomy_genus(taxonomy_dict):
    if GENUS in taxonomy_dict:
        genus = taxonomy_dict[GENUS]

    else:
        genus = 'NA'

    return genus


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
    output_file_name = INPUT_FILE_PREFIX + month_year_stamp + OUTPUT_FILE_SUFFIX + FILE_EXTENSION

    with open(output_file_name, 'w') as output_file:
        output_file.write('\t'.join([str(x) for x in ('Species', 'Num_assemblies', 'Kingdom', 'Phylum', 'Class', \
            'Order', 'Family', 'Genus')]) + '\n')

        LOG_FILE.write('\n#########################################################\n')
        LOG_FILE.write("Querying Species' taxonomy\n")
        LOG_FILE.write('#########################################################\n')

        num_species_retrieved = 0
        num_species_irretrievable = 0

        for i in range(len(species_list)):
            species_corrected = re.sub(r'(\[|\]|\(.+?\))', '', species_list[i])
            taxonomy_record = get_taxonomy_record(email, api_key, species_corrected)
            taxonomy_dict = get_full_taxonomy(taxonomy_record)
            kingdom = get_taxonomy_kingdom(taxonomy_dict)
            phylum = get_taxonomy_phylum(taxonomy_dict)
            taxonomy_class = get_taxonomy_class(taxonomy_dict)
            order = get_taxonomy_order(taxonomy_dict)
            family = get_taxonomy_family(taxonomy_dict)
            genus = get_taxonomy_genus(taxonomy_dict)

            if taxonomy_dict:
                num_species_retrieved += 1
                log_message = "Species " + str(i+1) + ": " + species_corrected + " taxonomy retrieved\n"
                print(log_message, end='')

            else:
                num_species_irretrievable += 1
                log_message = "Species " + str(i+1) + ": " + species_corrected + " taxonomy cannot be retrieved\n"
                print(log_message, end='')

            output_file.write('\t'.join([str(x) for x in (species_corrected ,num_assemblies[i], kingdom, \
                phylum, taxonomy_class, order, family, genus)]) + '\n')

    output_file.close()
    log_message_success = 'Total ' + str(num_species_retrieved) + " species' taxonomy retrieved successfully\n"
    LOG_FILE.write(log_message_success)
    log_message_failure = 'Total ' + str(num_species_irretrievable) + " species' taxonomy cannot be retrieved\n"
    LOG_FILE.write(log_message_failure)

if __name__ == '__main__':
	main()
