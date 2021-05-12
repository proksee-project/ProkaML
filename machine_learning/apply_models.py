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

    def __init__(self, dataframe, best_fit_model, training_taxon, testing_taxon):
        self.dataframe = dataframe
        self.best_fit_model = best_fit_model
        self.training_taxon = training_taxon
        self.testing_taxon = testing_taxon

    def apply_model_to_database(self):
        """
        Applies best fitting model to entire database

        PARAMETERS
            best_fit_model (obj.): object of class RandomForestClassifier with highest average cross-validation score

        RETURNS
            self.dataframe (obj.): object of class pandas.DataFrame having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 30 columns of different assembly attributes (str, int, float)
        """

        classifier = MachineLearningClassifier(self.dataframe, self.training_taxon)
        columns_of_interest = classifier.select_columns_of_interest(self.testing_taxon)
        prediction_column_label = 'train_' + self.training_taxon + '_test_' + self.testing_taxon + '_prob'
        database_dataframe = self.dataframe[columns_of_interest[:-1]]

        # Predicting assembly QC probability and assigning prediction probability column to database
        self.dataframe[prediction_column_label] = self.best_fit_model.predict_proba(database_dataframe)[:, 0]