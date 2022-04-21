import pandas as pd
import os
import re
from busco_lineage import busco_taxonomic_lists
import sys
filename = sys.argv[1]
output_file = open(filename.split('_')[0] + '_assembly_buscolineage.txt', 'w')

busco_kingdom, busco_phylum, busco_class, busco_order, busco_family = busco_taxonomic_lists()

contigs_list = os.listdir('predownloaded_assemblies')
file_dict = {}
for file in contigs_list:
	genbank_id_regex = re.search(r'(GCA_.+?)_.*', file)
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

input_file = open(filename, 'r')
for line in input_file:
	row = line.rstrip().split('\t')
	species = row[0]
	assembly_id = row[1]
	genbank_id = row[2]
	refseq_id = row[3]
	refseq_exclusion_reason = row[4]

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
				print('\t'+ "No busco lineage")
		else:
			print('\t'+ "taxonomy not found\n")
