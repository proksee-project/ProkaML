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
from pathlib import Path
from clean_metadata import OrganizeMetadata
from preprocess_metadata import SpeciesNormalization

START_DIR = Path(__file__).resolve().parents[1]
metadata_file = '{}/well_represented_species_metadata_added_attributes.txt'.format(str(START_DIR))

starting_dataframe = pd.read_csv(metadata_file, low_memory=False, sep='\t')
organize_dataframe = OrganizeMetadata(starting_dataframe)
prenormalized_dataframe = organize_dataframe.prenormalize_metadata()
print('Pre-normalization step complete. Assembly methods and sequencing technologies are curated. ', end='')
print('Long read assemblies have been excluded')

species_normalization = SpeciesNormalization(prenormalized_dataframe)
species_normalization.apply_normalization_to_database()
print('Post-normalization step complete. Metadata with appended normalized attributes ', end='')
print("are written to 'well_represented_species_metadata_normalized.txt'")
print("Database of species specific median attributes are written to 'species_median_log_metrics.txt'")
