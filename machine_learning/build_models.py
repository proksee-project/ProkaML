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
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Defining constants for RandomForestClassifier parameters
NUM_DECISION_TREES = 100
DECISION_CRITERION = "gini"
MAX_DEPTH = None # None translates to unlimited depth for seed decision trees
FEATURE_SELECTION_CRITERION = "auto"
RANDOM_SEED = 42

# Defining instance of RandomForestClassifier
rfc = RandomForestClassifier(n_estimators=NUM_DECISION_TREES,
                             criterion=DECISION_CRITERION,
                             max_depth=MAX_DEPTH,
                             max_features=FEATURE_SELECTION_CRITERION,
                             random_state=RANDOM_SEED)


class MachineLearningClassifier():
    """
    A class for generating and evaluating machine learning models from assembly QC data

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 29 columns of different assembly attributes (str, int, float)
    """

    def __init__(self, dataframe):
        """
        Initializes the MachineLearningClassifier class

        PARAMETERS
            dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 29 columns of different assembly attributes (str, int, float)
        """

        self.dataframe = dataframe

    def select_columns_of_interest(self, taxonomy_resolution):

        logn50_normalized = 'logn50_normalized_' + str(taxonomy_resolution) 
        logcontigcount_normalized = 'logcontigcount_normalized_' + str(taxonomy_resolution)
        logl50_normalized = 'logl50_normalized_' + str(taxonomy_resolution)
        logtotlen_normalized = 'logtotlen_normalized_' + str(taxonomy_resolution)
        gccontent_normalized = 'gccontent_normalized_' + str(taxonomy_resolution)
        label = 'label'

        columns_of_interest = [logn50_normalized, logcontigcount_normalized, logl50_normalized, logtotlen_normalized,
                               gccontent_normalized, label]

        return columns_of_interest

    def curate_labels(self):
        """
        Curates labels based on assembly inclusion or exclusion in NCBI RefSeq database

        RETURNS
           inclusion_data, exclusion_data : tuple of dataframes (obj.) for inclusion and exclusion datasets. Objects are of
           class pandas.Dataframe with two-dimensional data structures of ~95,000 rows and 6 columns (float, int) for inclusion
           dataset and ~18,000 rows and 6 columns (float, int) for exclusion dataset 
        """

        REFSEQ_ACCESSION = 'Refseq Accession'
        REFSEQ_EXCLUSION = 'Refseq Exclusion Reason'
        inclusion_label = 1
        exclusion_label = 2

        RefSeq_included_data = self.dataframe[self.dataframe[REFSEQ_ACCESSION].notnull()]
        RefSeq_included_data = RefSeq_included_data.assign(label=inclusion_label)

        RefSeq_excluded_data = self.dataframe[self.dataframe[REFSEQ_ACCESSION].isnull()]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_data[~RefSeq_excluded_data[REFSEQ_EXCLUSION].
                                                             str.contains(r'multi-isolate|\[\]', regex=True)]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_multi_isolate.assign(label=exclusion_label)

        # Creating dataframes of inclusion and exclusion datasets with columns of interest
        taxonomy_resolution = 'species'
        columns_of_interest = self.select_columns_of_interest(taxonomy_resolution)
        inclusion_data = RefSeq_included_data[columns_of_interest]
        exclusion_data = RefSeq_excluded_multi_isolate[columns_of_interest]

        return inclusion_data, exclusion_data

    def generate_models(self, inclusion_data, exclusion_data, num_iterations):
        """
        Generates random forests models and calculates cross-validation scores

        PARAMETERS
            inclusion_data (obj.): Object of class pandas.Dataframe with ~95,000 rows and 6 columns (float, int)
            exclusion_data (obj.): Object of class pandas.Dataframe with ~18,000 rows and 6 columns (float, int)
            num_iterations (int): number of random forests models to be generated

        RETURNS
            list_models_and_scores (list): list of sublists, with each sublist containing an object of class
            RandomForestClassifier and average cross-validation score (float)
        """

        list_models_and_scores = []
        N_FOLD_CV = 10
        taxonomy_resolution = 'species'

        for i in range(0, num_iterations):

            # Randomized balanced subsamping of inclusion data
            inclusion_data_balanced = inclusion_data.sample(n=exclusion_data.shape[0])
            integrated_data = pd.concat([inclusion_data_balanced, exclusion_data], ignore_index=True)

            # Defining feature space and label vector
            feature_space = integrated_data.drop(self.select_columns_of_interest(taxonomy_resolution)[-1], axis=1)
            label_vector = integrated_data[self.select_columns_of_interest(taxonomy_resolution)[-1]]

            # Fitting random forests model to the training data
            rfc.fit(feature_space, label_vector)
            print('Fitting random forest model {}. Evaluating cross-validation score'.format(int(i+1)))

            # Evaluating models by N (N=10) fold cross-validation
            model_average_cv_score = np.average(cross_val_score(rfc, feature_space, label_vector, cv=N_FOLD_CV))

            list_models_and_scores.append([rfc, model_average_cv_score])

        return list_models_and_scores

    def return_best_model(self):
        """
        Returns prediction model with the highest average cross-validation score

        RETURNS
            best_fit_model (obj.): object of class RandomForestClassifier with highest average cross-validation score
        """

        inclusion_data, exclusion_data = self.curate_labels()
        num_iterations = 10
        list_models_and_scores = self.generate_models(inclusion_data, exclusion_data, num_iterations)
        CLASSIFIER_MODEL_INDEX = 0

        max_score = max(list_models_and_scores[i][1] for i in range(0, len(list_models_and_scores)))
        max_index = [list_models_and_scores[i][1] for i in range(0, len(list_models_and_scores))].index(max_score)
        print('Random Forest model {} has highest average cross-validation score. '.format(int(max_index+1)), end='')
        print('Returning model {} as a python object.'.format(int(max_index+1)))
        best_fit_model = list_models_and_scores[max_index][CLASSIFIER_MODEL_INDEX]

        return best_fit_model

    def apply_model_to_database(self, best_fit_model, taxonomy_resolution):
        """
        Applies best fitting model to entire database

        PARAMETERS
            best_fit_model (obj.): object of class RandomForestClassifier with highest average cross-validation score

        RETURNS
            self.dataframe (obj.): object of class pandas.DataFrame having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 30 columns of different assembly attributes (str, int, float)
        """

        PREDICTION_COLUMN = 'RefSeq_'+ str(taxonomy_resolution) + '_prob'
        columns_of_interest = self.select_columns_of_interest(taxonomy_resolution)
        database_dataframe = self.dataframe[columns_of_interest[:-1]]

        # Predicting assembly QC probability and assigning prediction probability column to database
        self.dataframe[PREDICTION_COLUMN] = best_fit_model.predict_proba(database_dataframe)[:, 0]
