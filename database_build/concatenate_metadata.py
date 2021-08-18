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
from datetime import date
import glob
import re
import os

THRESHOLD_MAJOR = 1000
PREFIX_MAJOR = 'major'
THRESHOLD_LARGE = 100
PREFIX_LARGE = 'large'
THRESHOLD_INTERMEDIATE = 10
PREFIX_INTERMEDIATE = 'intermediate'

METADATA_DIRECTORY = 'species_metadata'
METADATA_FILE_EXTENSION = '_chunk.+?_metadata.txt'
JOIN_CHARACTER = '_'
SPLIT_CHARACTER = ' '
SEPARATOR = '\t'
METADATA_FILE_LIST = os.listdir(METADATA_DIRECTORY)
METADATA_COLUMNS = ['Organism Name', 'Strain/Isolate', 'Assembly Name', 'Genbank Accession',
					'Refseq Accession', 'Genome Coverage', 'Submission Date', 'Last Update Date',
					'Refseq Exclusion Reason', 'ContigN50', 'Contig count', 'ContigL50',
					'Total length', 'Assembly Method', 'Sequencing Technology']
SPECIES = 'Species'
KINGDOM = 'Kingdom'
NUM_ASSEMBLIES = 'Num_assemblies'
METADATA_FILES = 'Metadata_files'
INTEGRATED_METADATA_FILE = 'well_represented_species_metadata.txt'

def threshold_metadata(dataframe):
	dataframe = dataframe[dataframe[KINGDOM].str.contains(r'Archaea|Bacteria|NA\s\(Cannot')]
	df_major = dataframe[dataframe[NUM_ASSEMBLIES] >= THRESHOLD_MAJOR]
	df_large = dataframe[(dataframe[NUM_ASSEMBLIES] >= THRESHOLD_LARGE) & (dataframe[NUM_ASSEMBLIES] < THRESHOLD_MAJOR)]
	df_intermediate = dataframe[(dataframe[NUM_ASSEMBLIES] >= THRESHOLD_INTERMEDIATE) & (dataframe[NUM_ASSEMBLIES] < THRESHOLD_LARGE)]
	categorical_list_dataframes = [df_major, df_large, df_intermediate]

	return categorical_list_dataframes

def get_metadata_files(species):
	species_pattern = JOIN_CHARACTER.join(species.split(SPLIT_CHARACTER)) + METADATA_FILE_EXTENSION
	species_file_list = list(filter(re.compile(species_pattern).match, METADATA_FILE_LIST))
	species_file_paths_list = [os.path.join(METADATA_DIRECTORY, file) for file in species_file_list]

	return species_file_paths_list

def merge_series_list(metadata_file_series):
	"""
	Converts series (containing list of individual lists) to a single merged list
	"""
	merged_file_list = []
	for species_file_list in metadata_file_series:
		if len(species_file_list) > 0:
			for i in species_file_list:
				merged_file_list.append(i)

	return merged_file_list


def main():
	month_year_stamp = date.today().strftime("%b_%Y")
	assembly_count_dataframe = pd.read_csv('species_assemblycounts_' + month_year_stamp + '_taxonomy.txt', sep='\t')
	categorical_list_dataframes = threshold_metadata(assembly_count_dataframe)
	all_species_integrated_dataframe_list = []
	
	for i in range(0, len(categorical_list_dataframes)):
		if i == 0:
			category = PREFIX_MAJOR
		if i == 1:
			category = PREFIX_LARGE
		if i == 2:
			category = PREFIX_INTERMEDIATE

		categorical_integrated_metadata_file = category + '_species_metadata.txt'

		df = categorical_list_dataframes[i]
		df[METADATA_FILES] = df.apply(lambda row: get_metadata_files(row[SPECIES]), axis=1)
		categorical_metadata_file_list = merge_series_list(df[METADATA_FILES].to_list())
		categorical_metadataframe_list = []
		for j in categorical_metadata_file_list:
			species_chunk_dataframe = pd.read_csv(j, sep=SEPARATOR, names=METADATA_COLUMNS)
			categorical_metadataframe_list.append(species_chunk_dataframe)

		categorical_integrated_dataframe = pd.concat(categorical_metadataframe_list)
		categorical_integrated_dataframe.to_csv(categorical_integrated_metadata_file, sep=SEPARATOR, index=False)

		all_species_integrated_dataframe_list.append(categorical_integrated_dataframe)

	all_species_integrated_dataframe = pd.concat(all_species_integrated_dataframe_list)
	all_species_integrated_dataframe.to_csv(INTEGRATED_METADATA_FILE, sep=SEPARATOR, index=False)

if __name__ == '__main__':
	main()
