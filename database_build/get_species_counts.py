import pandas as pd
import constants as const
import os
from collections import defaultdict
import re
from get_species_taxonomy import Taxonomy
import argparse
from datetime import date

num_assemblies = defaultdict(int)
entrez_metadata = defaultdict(list)
filelist = os.listdir(const.ENTREZ_METADATA_DIR)

for file in filelist:
	filepath = os.path.join(const.ENTREZ_METADATA_DIR, file)
	df = pd.read_csv(filepath, sep=const.SEPARATOR, header=None, names=const.METADATA_COLUMNS)
	d = dict(tuple(df.groupby(['Organism Name'])))
	for species, dataframe in d.items():
		num_assemblies[species] += dataframe.shape[0]
		entrez_metadata[species].append(dataframe)

reorganized_metadata_dir = 'species_reorganized_metadata'
if not os.path.exists(reorganized_metadata_dir):
	os.mkdir(reorganized_metadata_dir)

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
output_file_name = const.TAXONOMY_FILE_PREFIX + month_year_stamp + const.TAXONOMY_FILE_SUFFIX + \
    const.FILE_EXTENSION
unexpected_chars = re.compile('[^\sa-zA-Z0-9\.-]')

output_file = open(output_file_name, const.WRITE_MODE)
output_file.write('\t'.join([str(x) for x in ('Species', 'Num_assemblies', 'Kingdom', 'Phylum', 'Class', \
	'Order', 'Family', 'Genus')]) + '\n')

for species, dataframe_list in entrez_metadata.items():
	df = pd.concat(entrez_metadata[species])
	species_corrected =	re.sub(unexpected_chars, const.EMPTY_STRING, species)

	if num_assemblies[species] >= const.ASSEMBLY_COUNT_LOWERBOUND:
		taxonomy = Taxonomy(species_corrected, email, api_key)
		taxonomy_list = taxonomy.get_full_taxonomy()
		output_file.write('\t'.join([species_corrected, str(num_assemblies[species]), *taxonomy_list]) + '\n')
		species_metadata_file = '_'.join(species_corrected.split(' ')) + const.METADATA_SUFFIX + const.FILE_EXTENSION
		df.to_csv(os.path.join(reorganized_metadata_dir, species_metadata_file), sep=const.SEPARATOR, index=False)
