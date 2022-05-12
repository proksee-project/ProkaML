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
import os
import re
from get_busco_lineage import busco_taxonomic_lists
import argparse


def main():

	my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] input_file',
                                    description='Maps genbank assemdly ID to assembly file')
	my_parser.add_argument('input_file',
                        type=str,
                        help='file containing assembly metadata')

	args = my_parser.parse_args()

	input_file = args.input_file
	output_file = open(input_file.split('.')[0] + '_buscolineage.' + input_file.split('.')[1], 'w')

	busco_kingdom, busco_phylum, busco_class, busco_order, busco_family = busco_taxonomic_lists()

	contigs_list = os.listdir('additional_species_metadata')
	file_dict = {}
	for file in contigs_list:
		genbank_id_regex = re.search(r'(GCA_.+?)_.*?\.gz$', file)
		if genbank_id_regex is not None:
			file_dict[genbank_id_regex.group(1)] = file 


	taxonomy_dict = {}
	with open('species_taxonomy.txt','r') as proka_taxonomy_db:
		next(proka_taxonomy_db)
		for line in proka_taxonomy_db:
			row = line.rstrip().split('\t')
			species = row[0]
			kingdom = row[2].lower()
			phylum = row[3].lower()
			taxonomy_class = row[4].lower()
			order = row[5].lower()
			family = row[6].lower()
			genus= row[7].lower()
			taxonomy_dict[species] = [family, order, taxonomy_class, phylum, kingdom]

	input = open(input_file, 'r')
	line_count = 0
	for line in input:
		row = line.rstrip().split('\t')
		
		line_count += 1
		if line_count == 1:
			output_file.write(line.rstrip() + '\tAssembly filename\tBusco Lineage\tLineage Type\n')

		species = row[0]
		genbank_id = row[3]

		if genbank_id in file_dict:
			output_file.write('\t'.join(row) + '\t' + file_dict[genbank_id] + '\t')
			if species in taxonomy_dict:
				if taxonomy_dict[species][0] in busco_family:
					output_file.write(taxonomy_dict[species][0] + '\t' + 'Family\n')
				elif taxonomy_dict[species][1] in busco_order:
					output_file.write(taxonomy_dict[species][1] + '\t' + 'Order\n')
				elif taxonomy_dict[species][2] in busco_class:
					output_file.write(taxonomy_dict[species][2] + '\t' + 'Class\n')
				elif taxonomy_dict[species][3] in busco_phylum:
					output_file.write(taxonomy_dict[species][3] + '\t' + 'Phylum\n')
				elif taxonomy_dict[species][4] in busco_kingdom:
					output_file.write(taxonomy_dict[species][4] + '\t' + 'Kingdom\n')
				else:
					output_file.write('\t'+ "No busco lineage")
			else:
				output_file.write('\t'+ "taxonomy not found\n")

if __name__ == '__main__':
	main()
