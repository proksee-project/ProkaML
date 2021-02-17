'''
Copyright:
University of Manitoba & National Microbiology Laboratory, Canada, 2020
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

'''
This program processes NCBI annotated metadata to calculate median of log of microbial specific
species specific genomic attributes
'''

import pandas as pd
import numpy as np
import os
import sys
import re 

#function to exclude species without valid taxonomy names
def exclude_species_regex():
	#list of regular expressions to exclude invalid species taxonomy names, can be expanded
	exclude_species_regex = ['uncultured', 
							'metagenome',
							'^bacterium$'
							]
	exclude_species_regex_compile = "|".join(exclude_species_regex)
	exclude_regexes = re.compile(exclude_species_regex_compile)

	return exclude_regexes

#function to read downloaded metadata tsv files
def species_to_process(filepath, exclude_regexes):
	df = pd.read_csv(filepath, sep='\t')

	#exclude species without valid taxonomy names and avoiding pandas append an 'Unnamed' column 
	species_ok = df[~df['Organism Name'].str.contains(exclude_regexes)]
	species_ok = species_ok.loc[:, ~species_ok.columns.str.contains('^Unnamed')]

	return species_ok

#function to group metadata by individual species
#takes integrated dataframe, splits by species and returns list of splitted dataframes per species 
def species_grouping(species_ok):
	species_grp = species_ok.groupby('Organism Name')
	species_grp_list = [species_grp.get_group(x) for x in species_grp.groups]

	return species_grp_list

#function to calculate median of log of numerical genomic attributes per species 
#attributes are N50, Contig count, L50, total length, coverage, can be expanded 
def median_log_species(df):
	df = df.assign(logn50=round(np.log10(df['ContigN50']),3))
	df = df.assign(logcontigcount=round(np.log10(df['Contig count']),3))
	df = df.assign(logl50=round(np.log10(df['ContigL50']),3))
	df = df.assign(logtotlen=round(np.log10(df['Total length']),3))
	df.loc[df['Genome Coverage'] > 0, 'adj coverage'] = df['Genome Coverage']
	df.loc[df['Genome Coverage'] == 0, 'adj coverage'] = 0.1
	df = df.assign(logcoverage=round(np.log10(df['adj coverage']),3))

	median_log_n50 = df['logn50'].median()
	median_log_contigcount = df['logcontigcount'].median()
	median_log_l50 = df['logl50'].median()
	median_log_totlen = df['logtotlen'].median()
	median_log_coverage = df['logcoverage'].median()

	return median_log_n50, median_log_contigcount, median_log_l50, median_log_totlen, median_log_coverage

#function to write median of log of metrics to file
def species_reference_median_log(indiv_sp_df, outfile):
	outfile.write('{}\t'.format(indiv_sp_df.iloc[0, 0]))
	for i in range(0, len(median_log_species(indiv_sp_df))):
		outfile.write('{}\t'.format(median_log_species(indiv_sp_df)[i]))
	outfile.write('\n')


def main():
	#passing integrated metadata file as first argument
	filename = sys.argv[1]
	filedir = '/mnt/c/Users/arnab/Documents/GitHub/proksee-database/metadata/annotated_metadata/Original'
	filepath = os.path.join(filedir, filename)
	
	#excluding species without valid taxonomy names
	exclude_species_patterns = exclude_species_regex()
	species_accept_df = species_to_process(filepath, exclude_species_patterns)

	#splitting the resultant metadata frame per species
	species_grp_list = species_grouping(species_accept_df)

	#Output file header
	outfile_ref = open('species_median_log_metrics.txt', 'w')
	outfile_ref.write('Species\tlogn50\tlogcontigcount\tlogl50\tlogtotlen\tlogcoverage\n')
	
	#Writing median of log of calculated metrics per species
	for indiv_sp_df in species_grp_list:
		species_reference_median_log(indiv_sp_df, outfile_ref)

if __name__ == '__main__':
	main()
