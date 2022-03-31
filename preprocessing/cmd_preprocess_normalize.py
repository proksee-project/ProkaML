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

import os
import pandas as pd
from pathlib import Path
from clean_metadata import CleanMetadata
from calculate_log_median import MedianNormalization
from normalize_species import SpeciesNormalization
from normalize_genus import GenusNormalization


START_DIR = Path(__file__).resolve().parents[1]
ATTRIBUTE_DATABASE_FILE = 'well_represented_species_metadata_added_attributes.txt'
SEPARATOR = '\t'
LOW_MEMORY = False
INDEX = False
INDEX_COL = False
MEDIAN_DATABASE_FILE = 'species_normalization_factor_database.txt'
NORMALIZED_DATABASE_FILE = 'well_represented_species_metadata_normalized.txt'

starting_dataframe = pd.read_csv(os.path.join(START_DIR, ATTRIBUTE_DATABASE_FILE), low_memory=LOW_MEMORY, sep=SEPARATOR)
organize_dataframe = CleanMetadata(starting_dataframe)
cleaned_dataframe = organize_dataframe.clean_metadata()
print('Pre-normalization step complete.')
print('.....Assembly methods and sequencing technologies are curated.')
print('.....Long read assemblies are excluded.')
print('.....Species with invalid taxonomy names are excluded.')
print(".....Genus' column is appended to database.")

species_normalization = SpeciesNormalization(cleaned_dataframe)
with open(os.path.join(START_DIR, MEDIAN_DATABASE_FILE), 'w') as median_db_write:
    median_db_write.write('Species\tlogn50\tlogcontigcount\tlogl50\tlogtotlen\tlogcoverage\tgccontent\n')
    species_normalized_dataframe = species_normalization.apply_species_normalization(median_db_write)

print("Database of species specific median attributes are written to {}.".format(MEDIAN_DATABASE_FILE))
print("Species' normalized attributes are calculated.")
species_normalized_dataframe.to_csv(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), sep=SEPARATOR, index=INDEX)

'''
species_median_dataframe = pd.read_csv(os.path.join(START_DIR, MEDIAN_DATABASE_FILE), sep=SEPARATOR, index_col=INDEX_COL)
genus_normalization = GenusNormalization(species_normalized_dataframe)
with open(os.path.join(START_DIR, MEDIAN_DATABASE_FILE), 'a') as median_db_append:
    genus_normalization.write_median_values(species_median_dataframe, median_db_append)

genus_normalized_dataframe = genus_normalization.apply_genus_normalization(species_normalized_dataframe)
genus_normalized_dataframe.to_csv(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), sep=SEPARATOR, index=INDEX)
print("Database of genus specific median attributes are appended to {}".format(MEDIAN_DATABASE_FILE))
print("Genus' normalized attributes are calculated.")
print("Species' and genus' normalized attributes for entire database are written to {}".format(NORMALIZED_DATABASE_FILE))
'''