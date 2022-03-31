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
from build_models import MachineLearningClassifier
from sklearn.ensemble import RandomForestClassifier


class MachineLearningPrediction():

    def __init__(self, model_build_dataframe, model_object, model_training_taxon, model_testing_taxon):
        self.model_build_dataframe = model_build_dataframe
        self.model_object = model_object
        self.model_training_taxon = model_training_taxon
        self.model_testing_taxon = model_testing_taxon

        self.classifier = MachineLearningClassifier(self.model_build_dataframe, self.model_training_taxon)
        self.columns_of_interest = self.classifier.select_columns_of_interest(self.model_testing_taxon)
        self.prediction_column_label = 'train_' + self.model_training_taxon + '_test_' + self.model_testing_taxon + '_prob'

    def apply_model(self, model_test_dataframe):

        test_feature_space = model_test_dataframe[self.columns_of_interest[:-1]]

        # Predicting assembly QC probability and assigning prediction probability column to database
        model_test_dataframe[self.prediction_column_label] = self.model_object.predict_proba(test_feature_space)[:, 0]
