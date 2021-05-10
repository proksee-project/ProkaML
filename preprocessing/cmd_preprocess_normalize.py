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
from clean_metadata import OrganizeMetadata
from preprocess_metadata import TaxonomicalNormalization

START_DIR = Path(__file__).resolve().parents[1]
metadata_file = '{}/well_represented_species_metadata_added_attributes.txt'.format(str(START_DIR))
NORMALIZED_DATABASE_FILE = open(os.path.join(START_DIR, 'norm_test.txt'), 'w')
SEPARATOR = '\t'


starting_dataframe = pd.read_csv(metadata_file, low_memory=False, sep=SEPARATOR)
organize_dataframe = OrganizeMetadata(starting_dataframe)
prenormalized_dataframe = organize_dataframe.prenormalize_metadata()
prenormalized_dataframe.to_csv(open(os.path.join(START_DIR, 'temp.txt'), 'w'), sep=SEPARATOR, mode='w', index=False)
print('Pre-normalization step complete. Assembly methods and sequencing technologies are curated. ', end='')
print('Long read assemblies have been excluded')

taxonomical_normalization = TaxonomicalNormalization(prenormalized_dataframe)
species_normalized_dataframe = taxonomical_normalization.apply_species_normalization_to_database()
species_normalized_dataframe.to_csv(open(os.path.join(START_DIR, 'temp1.txt'), 'w'), sep=SEPARATOR, mode='w', index=False)

genus_normalized_dataframe = taxonomical_normalization.apply_genus_normalization_to_database()
genus_normalized_dataframe.to_csv(open(os.path.join(START_DIR, 'temp2.txt'), 'w'), sep=SEPARATOR, mode='w', index=False)

BRIDGING_COLUMN = 'Genbank Accession'
MERGE_TYPE = 'left'
COLUMNS_TO_MERGE = [BRIDGING_COLUMN,'logn50_normalized_genus', 'logcontigcount_normalized_genus', 'logl50_normalized_genus',
                   'logtotlen_normalized_genus', 'logcoverage_normalized_genus', 'gccontent_normalized_genus']
merged_normalized_dataframe = species_normalized_dataframe.merge(genus_normalized_dataframe[COLUMNS_TO_MERGE], on=BRIDGING_COLUMN, how=MERGE_TYPE)
merged_normalized_dataframe.to_csv(NORMALIZED_DATABASE_FILE, sep=SEPARATOR, mode='w', index=False)

print('Post-normalization step complete. Metadata with appended normalized attributes ', end='')
print("are written to 'well_represented_species_metadata_normalized.txt'")
print("Database of species specific median attributes are written to 'species_median_log_metrics.txt'")
