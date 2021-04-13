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
import re
import os
from pathlib import Path 
from clean_metadata import OrganizeMetadata
from grouping_normalizing import GroupedNormalizing

START_DIR = Path(__file__).resolve().parents[1]
metadata_file = '{}/well_represented_species_metadata_added_attributes.txt'.format(str(START_DIR))
median_output_file = open('species_median_log_metrics.txt', 'w')

dataframe = pd.read_csv(metadata_file, sep='\t')
preprocessing = OrganizeMetadata()

dataframe_assembly_cleaned = preprocessing.organize_assembly_method(dataframe)
dataframe_sequencing_platform_cleaned = preprocessing.organize_sequencing_technology(dataframe_assembly_cleaned)
long_read_index = preprocessing.subset_long_read(dataframe_sequencing_platform_cleaned)
dataframe_sequencing_platform_cleaned.drop(long_read_index, inplace=True)

prenormalized_dataframe = dataframe_sequencing_platform_cleaned
grouped_normalization = GroupedNormalizing()
valid_species_list = grouped_normalization.filter_valid_species(prenormalized_dataframe)

median_output_file.write('Species\tlogn50\tlogcontigcount\tlogl50\tlogtotlen\tlogcoverage\tgccontent\n')

for individual_species_dataframe in valid_species_list:
	grouped_normalization.write_median_statistics(individual_species_dataframe, median_output_file)

