
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
	"""
	Groups assembly metadata by species name

	RETURNS:
		species_dataframe_dict (dict): dictionary of species (str) mapping to assembly metadata (2 dimensional dataframe
		with rows of assemblies and columns of genomic attributes)
	"""

	list_dataframes = []

	for file in os.listdir(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ENTREZ_METADATA_DIR)):
		filepath = os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ENTREZ_METADATA_DIR, file)
		df = pd.read_csv(filepath, sep=const.FileFormat.SEPARATOR, header=None, \
			names=const.Metadata.METADATA_COLUMN_HEADERS)
		list_dataframes.append(df)

	dataframes_concatenated = pd.concat(list_dataframes, ignore_index=True)
	species_dataframe_dict = dict(tuple(dataframes_concatenated.groupby([const.Metadata.METADATA_COLUMN_HEADERS\
		[const.Metadata.METADATA_INDEX_SPECIES]])))

	return species_dataframe_dict


def perform_species_examination(special_characters, species_dataframe_dict, taxonomy_output_file):
	"""
	Filters assemblies by lowerbound threshold and prokaryote taxonomy

	PARAMETERS:
		special_characters (regex): regular expression for special characters
		species_dataframe_dict (dict): dictionary of species (str) mapping to assembly metadata (2 dimensional dataframe
		with rows of assemblies and columns of genomic attributes)
		taxonomy_output_file (str): file to which full taxonomy lineage for species are written

	RETURNS:
		num_species_written (int): number of species written to output file
		num_assemblies_written (int): number of assemblies written to output file
		excluded_species (list): list of species and corresponding assembly metadata excluded from analysis
	"""

	num_assemblies_written = 0
	num_species_written = 0
	excluded_species = []
	for species, dataframe in species_dataframe_dict.items():
		species_corrected =	re.sub(special_characters, const.FileFormat.EMPTY_STRING, species)
		if check_assembly_counts(dataframe) and is_kingdom_prokaryote(species_corrected, dataframe, taxonomy_output_file):
			write_assembly_metadata(species_corrected, dataframe)
			num_species_written += 1
			num_assemblies_written += dataframe.shape[0]
		else:
			excluded_species.extend([species_corrected, dataframe])

	return num_species_written, num_assemblies_written, excluded_species


def check_assembly_counts(dataframe):
	"""
	Checks number of assemblies for a species

	PARAMETERS:
		dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
	
	RETURNS:
		(bool) whether number of assemblies are greater than a lowerbound threshold 
	"""

	if dataframe.shape[0] >= const.Assembly.ASSEMBLY_COUNT_LOWERBOUND:
		return True

	else:
		return False


def is_kingdom_prokaryote(species_corrected, dataframe, taxonomy_output_file):
	"""
	Checks whether species is a prokaryote

	PARAMETERS:
		species_corrected (str): species corrected for special characters
		dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
		taxonomy_output_file (str): file to which full taxonomy lineage for species are written
	
	RETURNS:
		(bool) whether species is a prokaryote
	"""
	
	my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key',
										description='Gets full taxonomical lineage for species')
	my_parser.add_argument('email', type=str, help='user email address')
	my_parser.add_argument('api_key', type=str, help='NCBI user API key')                      
	args = my_parser.parse_args()
	email = args.email
	api_key = args.api_key

	taxonomy_lineage = Taxonomy(species_corrected, email, api_key).fetch()
	if taxonomy_lineage[const.Taxonomy.SUPERKINGDOM_INDEX] in const.Taxonomy.PROKARYOTES:
		write_species_taxonomy(species_corrected, dataframe, taxonomy_lineage, taxonomy_output_file)
		return True

	else:
		return False


def write_species_taxonomy(species_corrected, dataframe, taxonomy_lineage, taxonomy_output_file):
	"""
	Writes full taxonomy lineage for a given species

	PARAMETERS:
		species_corrected (str): species corrected for special characters
		dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
		taxonomy_lineage (list): full taxonomy lineage for a given species
		taxonomy_output_file (str): file to which full taxonomy lineage for species are written
	
	POST:
		Writes species taxonomy lineage upon True evaluation of prokaryote
	"""

	taxonomy_output_file.write(const.FileFormat.SEPARATOR.join([species_corrected, str(dataframe.shape[0]), \
		*taxonomy_lineage]) + '\n')


def write_assembly_metadata(species_corrected, dataframe):
	"""
	Writes species' assembly metadata

	PARAMTERS:
		species_corrected (str): species corrected for special characters
		dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
	
	POST:
		Writes species' assembly metadata
	"""

	species_metadata_file = '_'.join(species_corrected.split(' ')) + const.Assembly.METADATA_SUFFIX + \
		const.FileFormat.TEXT
	dataframe.to_csv(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.REORGANIZED_METADATA_DIR, \
		species_metadata_file), sep=const.FileFormat.SEPARATOR, index=False)


def generate_taxonomy_file():
	"""
	Generates taxonomy output file based on current month and year

	RETURNS:
		taxonomy_output_file (str): the taxonomy output file
	"""

	month_year_stamp = date.today().strftime(const.FileFormat.DATE_FORMAT)
	taxonomy_output_file_name = const.Taxonomy.TAXONOMY_FILE_PREFIX + month_year_stamp + \
		const.Taxonomy.TAXONOMY_FILE_SUFFIX + const.FileFormat.TEXT
	taxonomy_output_file = open(taxonomy_output_file_name, const.FileFormat.WRITE_MODE)
	taxonomy_output_file.write(const.FileFormat.SEPARATOR.join(const.Taxonomy.TAXONOMY_WRITE_COLUMNS) + '\n')

	return taxonomy_output_file


def write_excluded_metata(excluded_species, excluded_metadata_files):
	"""
	Writes metadata for species excluded from further analysis

	PARAMETERS:
		excluded_species (list): list of species and corresponding assembly metadata excluded from analysis
		excluded_metadata_files (list): list containing filenames (str) for excluded species and assemblies

	RETURNS:
		num_excluded_species (int): number of species excluded from further analysis
		num_excluded_assemblies (int): number of assemblies excluded from further analysis
	"""

	num_excluded_species = 0
	num_excluded_assemblies = 0
	with open(os.path.join(const.FileDirectories.DATABASE_PATH, excluded_metadata_files[0]), const.FileFormat.WRITE_MODE) as exclude:
		exclude.write(const.Assembly.COL_SPECIES + const.FileFormat.SEPARATOR + const.Assembly.COL_NUM_ASSEMBLIES + '\n')
		for i in range(0, len(excluded_species), 2):
			exclude.write(excluded_species[i] + const.FileFormat.SEPARATOR + str(excluded_species[i + 1].shape[0]) + '\n')
			if i == 0:
				mode = const.FileFormat.WRITE_MODE
				header = True
			else:
				mode = const.FileFormat.APPEND_MODE
				header = False
			excluded_species[i + 1].to_csv( os.path.join(const.FileDirectories.DATABASE_PATH, excluded_metadata_files[1]), mode=mode, \
				header=header, sep=const.FileFormat.SEPARATOR, index=False)
			num_excluded_species += 1
			num_excluded_assemblies += excluded_species[i + 1].shape[0]

	return num_excluded_species, num_excluded_assemblies


def main():
	"""
	Groups assembly metadata by species, examines assembly counts and taxonomy, generates report
	"""

	if not os.path.exists(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.REORGANIZED_METADATA_DIR)):
		os.mkdir(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.REORGANIZED_METADATA_DIR))
	if not os.path.exists(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ADDITIONAL_METADATA_DIR)):
		os.mkdir(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ADDITIONAL_METADATA_DIR))

	species_dataframe_dict = group_metadata()
	taxonomy_output_file = generate_taxonomy_file()
	special_characters = re.compile('[^\sa-zA-Z0-9\.-]')

	num_species_written, num_assemblies_written, excluded_species = perform_species_examination(special_characters,\
		species_dataframe_dict, taxonomy_output_file)
	excluded_metadata_files = [const.Assembly.EXCLUDED_SPECIES + const.FileFormat.TEXT, \
		const.Assembly.EXCLUDED_ASSEMBLIES + const.Assembly.METADATA_SUFFIX + const.FileFormat.TEXT]

	num_excluded_species, num_excluded_assemblies = write_excluded_metata(excluded_species, excluded_metadata_files)

	log = open(const.LogFiles.MAIN_LOG + const.FileFormat.TEXT, mode=const.FileFormat.APPEND_MODE)
	log.write('\n#########################################################\n')
	log.write('Reorganizing NCBI metadata by species, filtering by assembly counts and prokaryote taxonomy\n')
	log.write('#########################################################\n')
	log.write('Total {} assembly metadata for {} species written in species-specific files in format {}/'.format(\
		num_assemblies_written, num_species_written, os.path.join(os.path.basename(const.FileDirectories.DATABASE_PATH), \
		const.FileDirectories.REORGANIZED_METADATA_DIR)) + '[species]' + const.Assembly.METADATA_SUFFIX + const.FileFormat.TEXT + '\n')
	log.write('Total {} assembly metadata for {} species excluded from further analysis.\n'.format(num_excluded_assemblies,\
		num_excluded_species) + 'Excluded species with assembly counts written to {}\n'.format(os.path.join(os.path.basename(\
		const.FileDirectories.DATABASE_PATH), excluded_metadata_files[0])) + 'Excluded assembly metadata written to {}\n'.format(\
		os.path.join(os.path.basename(const.FileDirectories.DATABASE_PATH), excluded_metadata_files[1])))
	log.close()


if __name__ == '__main__':
	main()
