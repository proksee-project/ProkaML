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
SEPARATOR = '\t'
NORMALIZED_DATABASE_FILE = 'well_represented_species_metadata_normalized.txt'
DATABASE_SIGNATURE = '_well_represented'
model_build_dataframe = pd.read_table(os.path.join(START_DIR, NORMALIZED_DATABASE_FILE), low_memory=False)

print("Building random forest model for species' normalized attributes")
training_taxon = 'species'
testing_taxon = 'species'
classifier_species = MachineLearningClassifier(model_build_dataframe, training_taxon)
model_object = classifier_species.return_best_model(START_DIR, DATABASE_SIGNATURE)
print("Predicting on species' normalized attributes")
prediction_species = MachineLearningPrediction(model_build_dataframe, model_object, training_taxon, testing_taxon)
prediction_species.apply_model(model_test_dataframe)
