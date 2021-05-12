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
from pathlib import Path
from build_models import MachineLearningClassifier
from apply_models import MachineLearningPrediction
import joblib

START_DIR = Path(__file__).resolve().parents[1]
NORMALIZED_DATABASE_FILE = 'well_represented_species_metadata_normalized.txt'
#MODEL_FILENAME = 'random_forest_n50_contigcount_l50_totlen_gccontent.joblib.gz'

SPECIES_TAXON = 'species'
GENUS_TAXON = 'genus'
MODEL_FILENAME_EXTENSION = '_assemblyQC.joblib.gz'
COMPRESSION_TYPE = 'gzip'
COMPRESSION_LEVEL = 3
#PROBABILITY_DATABASE_FILE = 'well_represented_species_metadata_qc_probabilities.txt'
PROBABILITY_DATABASE_FILE = 'temp_prob.txt'
SEPARATOR = '\t'

filepath = str('{}/' + NORMALIZED_DATABASE_FILE).format(str(START_DIR))
dataframe = pd.read_table(filepath, low_memory=False)

print('Building random forest model for species normalized attributes')
classifier_species = MachineLearningClassifier(dataframe, SPECIES_TAXON)
species_best_fit_model = classifier_species.return_best_model()
model_filename = str(SPECIES_TAXON) + str(MODEL_FILENAME_EXTENSION)

# Exporting best fitting model as a *.joblib python object
joblib.dump(species_best_fit_model,
            os.path.join(START_DIR, model_filename),
            compress=(COMPRESSION_TYPE, COMPRESSION_LEVEL))

print('Building random forest model for genus normalized attributes')
classifier_genus = MachineLearningClassifier(dataframe, GENUS_TAXON)
genus_best_fit_model = classifier_genus.return_best_model()
model_filename = str(GENUS_TAXON) + str(MODEL_FILENAME_EXTENSION)

# Exporting best fitting model as a *.joblib python object
joblib.dump(genus_best_fit_model,
            os.path.join(START_DIR, model_filename),
            compress=(COMPRESSION_TYPE, COMPRESSION_LEVEL))

training_taxon1 = SPECIES_TAXON
testing_taxon1 = SPECIES_TAXON
prediction1 = MachineLearningPrediction(dataframe, species_best_fit_model, training_taxon1, testing_taxon1)
prediction1.apply_model_to_database()

training_taxon2 = SPECIES_TAXON
testing_taxon2 = GENUS_TAXON
prediction2 = MachineLearningPrediction(dataframe, species_best_fit_model, training_taxon2, testing_taxon2)
prediction2.apply_model_to_database()

training_taxon3 = GENUS_TAXON
testing_taxon3 = GENUS_TAXON
prediction3 = MachineLearningPrediction(dataframe, species_best_fit_model, training_taxon3, testing_taxon3)
prediction3.apply_model_to_database()

training_taxon4 = GENUS_TAXON
testing_taxon4 = SPECIES_TAXON
prediction4 = MachineLearningPrediction(dataframe, species_best_fit_model, training_taxon4, testing_taxon4)
prediction4.apply_model_to_database()

dataframe.to_csv(os.path.join(START_DIR, PROBABILITY_DATABASE_FILE),
                              sep=SEPARATOR, mode='w', index=False)
print('Assembly QC prediction probabilities appended to database')
