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

from gc_content import GCContentCalculate
import pandas as pd
import numpy as np
import subprocess
import gzip
import os
import glob
import shutil
import re
import argparse


def get_genus(genus_dict, species):
	if species in genus_dict and genus_dict[species] is not np.nan:
		genus = genus_dict[species]
	else:
		genus = np.nan

	return genus


def unzip_file(zipped_file_path):

	# Stripping .gz suffix from file name
	unzipped_file_path = zipped_file_path[:len(zipped_file_path) - 3]
	with gzip.open(zipped_file_path, 'rt') as zipped:
		with open(unzipped_file_path, 'w') as unzipped:
			shutil.copyfileobj(zipped, unzipped)
	
	return unzipped_file_path


def get_num_CDS(prokka_stats):
	regex_cds = re.search(r'CDS:\s(\d+)',prokka_stats)
	if regex_cds is not None:
		num_cds = regex_cds.group(1)
	else:
		num_cds = float('NaN')

	return num_cds


def get_num_tRNA(prokka_stats):
	regex_trna = re.search(r'tRNA:\s(\d+)',prokka_stats)
	if regex_trna is not None:
		num_trna = regex_trna.group(1)
	else:
		num_trna = float('NaN')

	return num_trna


def main():
	species_taxonomy_df = pd.read_csv('species_taxonomy.txt', sep='\t')
	genus_dict = dict(zip(species_taxonomy_df['Species'], species_taxonomy_df['Genus']))

	my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key input_file',
                                    description='Downloads assembly from NCBI and runs prokka')
	my_parser.add_argument('email',
                        type=str,
                        help='user email address')
	my_parser.add_argument('api_key',
                        type=str,
                        help='NCBI user API key')
	my_parser.add_argument('input_file',
                        type=str,
                        help='file containing assembly metadata')

	args = my_parser.parse_args()

	email = args.email
	api_key = args.api_key
	input_file = args.input_file

	df = pd.read_csv(input_file, sep='\t')
	gc = GCContentCalculate(email, api_key, df)
	gb_id_list = df['Genbank Accession'].to_list()
	species_list = df['Organism Name'].to_list()
	num_CDS_list = []
	num_tRNA_list = []

	for i in range(len(gb_id_list)):
		assembly_file_path_zipped = gc.download_assembly(gb_id_list[i])
		assembly_file_path_unzipped = unzip_file(assembly_file_path_zipped)
		genus = get_genus(genus_dict, species_list[i])
		if genus is not np.nan:
			prokka_cmd = 'prokka ' + assembly_file_path_unzipped + ' --force --quiet --outdir ' + gb_id_list[i] + ' --genus ' + genus
		else:
			prokka_cmd = 'prokka ' + assembly_file_path_unzipped + ' --force --quiet --outdir ' + gb_id_list[i]

		subprocess.call(prokka_cmd, shell=True)

		prokka_output_file = glob.glob(os.path.join(gb_id_list[i],'*.txt'))[0]
		with open(prokka_output_file, 'r') as prokka_output:
			prokka_stats = prokka_output.read()

		num_CDS = get_num_CDS(prokka_stats)
		num_CDS_list.append(num_CDS)
		num_tRNA = get_num_tRNA(prokka_stats)
		num_tRNA_list.append(num_tRNA)

		shutil.rmtree(gb_id_list[i])

	df['num_CDS'] = pd.Series(num_CDS_list)
	df['num_tRNA'] = pd.Series(num_tRNA_list)

	prokka_output_file = input_file.split('.')[0] + '_prokka.' + input_file.split('.')[1]
	df.to_csv(prokka_output_file, sep='\t', index=False)

if __name__ == '__main__':
	main()
