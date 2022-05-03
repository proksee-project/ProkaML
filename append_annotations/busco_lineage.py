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
import re

species_taxonomy_df = pd.read_csv('species_taxonomy.txt', sep='\t')
#Needs to be mapped correctly to get_species_counts.py output from Assembly pipeline

kingdom = species_taxonomy_df['Kingdom'].dropna().str.lower().unique().tolist()
phylum = species_taxonomy_df['Phylum'].dropna().str.lower().unique().tolist()
taxonomy_class = species_taxonomy_df['Class'].dropna().str.lower().unique().tolist()
order = species_taxonomy_df['Order'].dropna().str.lower().unique().tolist()
family = species_taxonomy_df['Family'].dropna().str.lower().unique().tolist()

busco_lineage = open('lineages_list.2019-11-27.txt','r')
busco_family = []
busco_order = []
busco_class = []
busco_phylum = []
busco_kingdom = []

def busco_taxonomic_lists():
	for line in busco_lineage:
		root_lineage_regex =  re.search(r'^\s(\w.*?)_odb', line.rstrip())

		if root_lineage_regex is not None:
			root_lineage = root_lineage_regex.group(1)
			if root_lineage in kingdom:
				busco_kingdom.append(root_lineage)

		else:
			other_lineage_regex = re.search(r'^\s{2,50}-\s(.+?)_odb', line.rstrip())
			if other_lineage_regex is not None:
				other_lineage = other_lineage_regex.group(1)

				if other_lineage in phylum:
					busco_phylum.append(other_lineage)
				elif other_lineage in taxonomy_class:
					busco_class.append(other_lineage)
				elif other_lineage in order:
					busco_order.append(other_lineage)
				elif other_lineage in family:
					busco_family.append(other_lineage)

	return busco_kingdom, busco_phylum, busco_class, busco_order, busco_family
