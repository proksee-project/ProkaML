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

import pandas as pd
import numpy as np
import os
import sys
import re 

class GroupedNormalizing():

	def __init__(self):
		pass

	def filter_valid_species(self, dataframe):
		#list of regular expressions to exclude invalid species taxonomy names, can be expanded
		exclude_species_patterns = ['uncultured', 
									'metagenome',
									'^bacterium$'
								   ]
		exclude_species = re.compile("|".join(exclude_species_patterns))
		valid_species_dataframe = dataframe[~dataframe['Organism Name'].str.contains(exclude_species)]
		valid_species_group = valid_species_dataframe.groupby('Organism Name')
		valid_species_dataframe_list = [valid_species_group.get_group(x) for x in valid_species_group.groups]

		return valid_species_dataframe_list

	def calculate_median_log(self, dataframe):
		dataframe = dataframe.assign(logn50=round(np.log10(dataframe['ContigN50']),3))
		dataframe = dataframe.assign(logcontigcount=round(np.log10(dataframe['Contig count']),3))
		dataframe = dataframe.assign(logl50=round(np.log10(dataframe['ContigL50']),3))
		dataframe = dataframe.assign(logtotlen=round(np.log10(dataframe['Total length']),3))
		dataframe.loc[dataframe['Genome Coverage'] > 0, 'adj coverage'] = dataframe['Genome Coverage']
		dataframe.loc[dataframe['Genome Coverage'] == 0, 'adj coverage'] = 0.1
		dataframe = dataframe.assign(logcoverage=round(np.log10(dataframe['adj coverage']),3))

		median_log_n50 = dataframe['logn50'].median()
		median_log_contigcount = dataframe['logcontigcount'].median()
		median_log_l50 = dataframe['logl50'].median()
		median_log_totlen = dataframe['logtotlen'].median()
		median_log_coverage = dataframe['logcoverage'].median()
		median_gc_content = dataframe['GCcontent'].median()

		return median_log_n50, median_log_contigcount, median_log_l50, median_log_totlen, median_log_coverage, median_gc_content

	def write_median_statistics(self, dataframe, outfile):
		outfile.write('{}\t'.format(dataframe.iloc[0, 0]))
		for i in range(0, len(self.calculate_median_log(dataframe))):
			outfile.write('{}\t'.format(self.calculate_median_log(dataframe)[i]))
		outfile.write('\n')
