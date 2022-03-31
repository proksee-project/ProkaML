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
#MODEL_TEST_FILE = 'intermediate_species_metadata_conditional_normalized.txt'
#model_test_dataframe = pd.read_table(os.path.join(START_DIR, MODEL_TEST_FILE), low_memory=False)
#PROBABILITY_DATABASE_FILE = 'temp2_prob.txt'
SEPARATOR = '\t'

# Case 1 Train on major+large+intermediate species' attr, test on species' attr
NORMALIZED_DATABASE_FILE = 'well_represented_species_metadata_normalized.txt'
DATABASE_SIGNATURE = '_well_represented'
model_build_dataframe = pd.read_table(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), low_memory=False)

print("Building random forest model for species' normalized attributes from major+large+intermediate")
training_taxon = 'species'
testing_taxon = 'species'
classifier_species = MachineLearningClassifier(model_build_dataframe, training_taxon)
model_object = classifier_species.return_best_model(START_DIR, DATABASE_SIGNATURE)
print("Predicting on species' normalized attributes")
prediction_species = MachineLearningPrediction(model_build_dataframe, model_object, training_taxon, testing_taxon)
#prediction_species.apply_model(model_test_dataframe)

'''
# Case 2 Train on major+large+intermediate genus' attr, test on genus' attr, consider species present
NORMALIZED_DATABASE_FILE = 'well_represented_species_metadata_normalized.txt'
DATABASE_SIGNATURE = '_well_represented'
model_build_dataframe = pd.read_table(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), low_memory=False)

print("Building random forest model for genus' normalized attributes from major+large+intermediate")
training_taxon = 'genus'
testing_taxon = 'genus'
classifier_species = MachineLearningClassifier(model_build_dataframe, training_taxon)
model_object = classifier_species.return_best_model(START_DIR, DATABASE_SIGNATURE)
print("Predicting on genus' normalized attributes")
prediction_genus = MachineLearningPrediction(model_build_dataframe, model_object, training_taxon, testing_taxon)
prediction_genus.apply_model(model_test_dataframe)


# Case 3 Train on major+large species' attr, test on intermediate genus' attr, consider species absent
NORMALIZED_DATABASE_FILE = 'major_large_species_metadata_normalized.txt'
DATABASE_SIGNATURE = '_major_large'
model_build_dataframe = pd.read_table(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), low_memory=False)

print("Building random forest model for species' normalized attributes from major+large")
training_taxon = 'species'
testing_taxon = 'genus'
classifier_species = MachineLearningClassifier(model_build_dataframe, training_taxon)
model_object = classifier_species.return_best_model(START_DIR, DATABASE_SIGNATURE)
print("Predicting on genus' normalized attributes")
prediction_genus = MachineLearningPrediction(model_build_dataframe, model_object, training_taxon, testing_taxon)
prediction_genus.apply_model(model_test_dataframe)


# Case 4 Train on major+large genus' attr, test on intermediate genus' attr, consider species absent
NORMALIZED_DATABASE_FILE = 'major_large_species_metadata_normalized.txt'
DATABASE_SIGNATURE = '_major_large'
model_build_dataframe = pd.read_table(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), low_memory=False)

print("Building random forest model for genus' normalized attributes from major+large")
training_taxon = 'genus'
testing_taxon = 'genus'
classifier_species = MachineLearningClassifier(model_build_dataframe, training_taxon)
model_object = classifier_species.return_best_model(START_DIR, DATABASE_SIGNATURE)
print("Predicting on genus' normalized attributes")
prediction_genus = MachineLearningPrediction(model_build_dataframe, model_object, training_taxon, testing_taxon)
prediction_genus.apply_model(model_test_dataframe)
'''

#model_test_dataframe.to_csv(os.path.join(START_DIR, PROBABILITY_DATABASE_FILE),
#                              sep=SEPARATOR, mode='w', index=False)
