
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
import constants as const
import os
from collections import defaultdict
import re
from get_species_taxonomy import Taxonomy
import argparse
from datetime import date


def group_metadata():
	filelist = os.listdir(const.ENTREZ_METADATA_DIR)
	list_dataframes = []

	for file in filelist:
		filepath = os.path.join(const.ENTREZ_METADATA_DIR, file)
		df = pd.read_csv(filepath, sep=const.SEPARATOR, header=None, names=const.METADATA_COLUMNS)
		list_dataframes.append(df)

	dataframes_concat = pd.concat(list_dataframes, ignore_index=True)
	species_dataframe_dict = dict(tuple(dataframes_concat.groupby([const.METADATA_COLUMNS[0]])))

	return species_dataframe_dict


def iterate_species(species_dataframe_dict, taxonomy_output_file, email, api_key):
	excluded_species = []
	excluded_assemblies = []
	unexpected_chars = re.compile('[^\sa-zA-Z0-9\.-]')

	for species, dataframe in species_dataframe_dict.items():
		species_corrected =	re.sub(unexpected_chars, const.EMPTY_STRING, species)
		included, excluded = filter_assembly_count_taxonomy(species_corrected, dataframe, 
			excluded_species, excluded_assemblies, taxonomy_output_file, email, api_key)

		num_assemblies_include += included
		num_assemblies_exclude += excluded

	return num_assemblies_include, num_assemblies_exclude, excluded_species, excluded_assemblies


def filter_assembly_count_taxonomy(species_corrected, dataframe, excluded_species, excluded_assemblies, 
	taxonomy_output_file, email, api_key):

	included = 0
	excluded = 0

	if dataframe.shape[0] >= const.ASSEMBLY_COUNT_LOWERBOUND:
		check_prokaryote_taxonomy(species_corrected, taxonomy_output_file, dataframe, 
			included, excluded, excluded_species, excluded_assemblies, email, api_key)

	else:
		excluded = dataframe.shape[0]
		excluded_species.append(species_corrected)
		excluded_assemblies.append(dataframe[const.METADATA_COLUMNS[3]].to_list())

	return included, excluded


def check_prokaryote_taxonomy(species_corrected, taxonomy_output_file, dataframe, included, 
	excluded, excluded_species, excluded_assemblies, email, api_key):

	taxonomy = Taxonomy(species_corrected, email, api_key)
	taxonomy_list = taxonomy.get_full_taxonomy()

	if taxonomy_list[0] in const.PROKARYOTES:
		taxonomy_output_file.write('\t'.join([species_corrected, str(dataframe.shape[0]), *taxonomy_list]) + '\n')
		included = dataframe.shape[0]

		species_metadata_file = '_'.join(species_corrected.split(' ')) + const.METADATA_SUFFIX + const.FILE_EXTENSION
		dataframe.to_csv(os.path.join(const.REORGANIZED_METADATA_DIR, species_metadata_file), sep=const.SEPARATOR, index=False)

	else:
		excluded = dataframe.shape[0]
		excluded_species.append(species_corrected)
		excluded_assemblies.append(dataframe[const.METADATA_COLUMNS[3]].to_list())

	return included, excluded

def main():
	if not os.path.exists(const.REORGANIZED_METADATA_DIR):
		os.mkdir(const.REORGANIZED_METADATA_DIR)
	if not os.path.exists(const.ADDITIONAL_METADATA_DIR):
		os.mkdir(const.ADDITIONAL_METADATA_DIR)

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

	month_year_stamp = date.today().strftime(const.DATE_FORMAT)
	taxonomy_output_file_name = const.TAXONOMY_FILE_PREFIX + month_year_stamp + const.TAXONOMY_FILE_SUFFIX + \
    	const.FILE_EXTENSION
	taxonomy_output_file = open(taxonomy_output_file_name, const.WRITE_MODE)
	taxonomy_output_file.write('\t'.join([str(x) for x in ('Species', 'Num_assemblies', 'Kingdom', 'Phylum', 'Class', \
		'Order', 'Family', 'Genus')]) + '\n')

	species_dataframe_dict = group_metadata()
	
	num_assemblies_included_init = 0
	num_assemblies_excluded_init = 0

	
	num_assemblies_included, num_assemblies_excluded, excluded_species, excluded_assemblies = \
		iterate_species(species_dataframe_dict, taxonomy_output_file, email, api_key, num_assemblies_included, num_assemblies_excluded,
		excluded_species, excluded_assemblies)

	print(num_assemblies_included)
	print(num_assemblies_excluded)

	excl_file = open('excl.txt','w')
	for species in excluded_species:
		excl_file.write(species + '\n')
	for assembly_list in excluded_assemblies:
		for assembly in assembly_list:
			excl_file.write(assembly + '\n')

	'''
	log = open(const.LOG_FILE + const.FILE_EXTENSION, mode=const.APPEND_MODE)
	log.write('\n#########################################################\n')
	log.write("Reorganizing NCBI metadata, filtering by assembly counts, relevant species' taxonomy\n")
	log.write('#########################################################\n')
	log.write('Assemblies included: {}\n'.format(assemblies_included))
	log.write('Assemblies excluded: {}\n'.format(assemblies_excluded))
	'''


if __name__ == '__main__':
	main()
