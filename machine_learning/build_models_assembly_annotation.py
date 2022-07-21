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
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
import joblib


class MachineLearningClassifier():
    """
    A class for generating and evaluating machine learning models from assembly QC data

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 29 columns of different assembly attributes (str, int, float)
    """

    # Defining constant parameters for Random Forest classifier model training and validation
    NUM_DECISION_TREES = 100
    DECISION_CRITERION = "gini"
    MAX_DEPTH = None # None translates to unlimited depth for seed decision trees
    FEATURE_SELECTION_CRITERION = "auto"
    RANDOM_SEED = 42
    NUM_ITERATIONS = 10
    N_FOLD_CV = 10
    NUM_PARALLEL_JOBS = -1 # Should be set to number of desired parallel processors. -1 engages all processors
    ROUNDING_DIGITS = 2

    def __init__(self, dataframe, training_taxon):
        """
        Initializes the MachineLearningClassifier class

        PARAMETERS
            dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 29 columns of different assembly attributes (str, int, float)
        """

        self.dataframe = dataframe
        self.training_taxon = training_taxon

        # Defining instance of RandomForestClassifier
        self.rfc = RandomForestClassifier(n_estimators = self.NUM_DECISION_TREES,
                                          criterion = self.DECISION_CRITERION,
                                          max_depth = self.MAX_DEPTH,
                                          max_features = self.FEATURE_SELECTION_CRITERION,
                                          random_state = self.RANDOM_SEED)

    def __select_columns_of_interest(self, input_taxon):

        logn50_normalized = 'logn50_normalized_' + str(input_taxon) 
        logcontigcount_normalized = 'logcontigcount_normalized_' + str(input_taxon)
        logl50_normalized = 'logl50_normalized_' + str(input_taxon)
        logtotlen_normalized = 'logtotlen_normalized_' + str(input_taxon)
        gccontent_normalized = 'gccontent_normalized_' + str(input_taxon)

        label = 'label'

        columns_of_interest = [logn50_normalized, logcontigcount_normalized, logl50_normalized, logtotlen_normalized,
                               gccontent_normalized, label]

        return columns_of_interest

    def __curate_labels(self):
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

        RefSeq_included_data = self.dataframe[self.dataframe[REFSEQ_ACCESSION].notnull() & 
                                              self.dataframe[REFSEQ_EXCLUSION].str.contains(r'\[\]', regex=True)]
        RefSeq_included_data = RefSeq_included_data.assign(label=inclusion_label)

        RefSeq_excluded_data = self.dataframe[self.dataframe[REFSEQ_ACCESSION].isnull()]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_data[~RefSeq_excluded_data[REFSEQ_EXCLUSION].
                                                             str.contains(r'multi-isolate|\[\]', regex=True)]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_multi_isolate.assign(label=exclusion_label)

        # Creating dataframes of inclusion and exclusion datasets with columns of interest
        columns_of_interest = self.__select_columns_of_interest(self.training_taxon)
        inclusion_data = RefSeq_included_data[columns_of_interest]
        exclusion_data = RefSeq_excluded_multi_isolate[columns_of_interest]

        return inclusion_data, exclusion_data

    def __generate_models(self, inclusion_data, exclusion_data):
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

        for i in range(0, self.NUM_ITERATIONS):

            # Randomized balanced subsamping of inclusion data
            inclusion_data_balanced = inclusion_data.sample(n=exclusion_data.shape[0])
            integrated_data = pd.concat([inclusion_data_balanced, exclusion_data], ignore_index=True)

            # Defining feature space and label vector
            feature_space = integrated_data.drop(self.__select_columns_of_interest(self.training_taxon)[-1], axis=1)
            label_vector = integrated_data[self.__select_columns_of_interest(self.training_taxon)[-1]]

            # Fitting random forests model to the training data
            self.rfc.fit(feature_space, label_vector)
            #print('Fitting random forest model {} by {} fold cross validation. '.format(int(i+1), self.N_FOLD_CV), end='')

            # Evaluating models by N fold cross-validation
            random_forest_cross_validation_predict = cross_val_predict(self.rfc,
                                                                       feature_space,
                                                                       label_vector,
                                                                       cv=self.N_FOLD_CV,
                                                                       n_jobs=self.NUM_PARALLEL_JOBS)
            conf_mat = confusion_matrix(label_vector, random_forest_cross_validation_predict)
            performance_metrics = self.__get_performance_metrics(conf_mat)
            print('{}\t{}\t{}'.format(*performance_metrics))

            list_models_and_scores.append([self.rfc, performance_metrics])

        return list_models_and_scores

    def __get_performance_metrics(self, confusion_matrix):

        PERCENT_MULTIPLIER = 100
        true_positive = int(confusion_matrix[0, 0])
        true_negative = int(confusion_matrix[1, 1])
        false_negative = int(confusion_matrix[0, 1])
        false_positive = int(confusion_matrix[1, 0])

        sensitivity = round(float(true_positive/(true_positive + false_negative))*PERCENT_MULTIPLIER, self.ROUNDING_DIGITS)
        specificity = round(float(true_negative/(true_negative + false_positive))*PERCENT_MULTIPLIER, self.ROUNDING_DIGITS)
        average_sensitivity_specificy = round((sensitivity + specificity)/2, self.ROUNDING_DIGITS)
        accuracy = round(float((true_positive + true_negative)/
                               (true_positive + true_negative + false_positive + false_negative)
                              )*PERCENT_MULTIPLIER, self.ROUNDING_DIGITS
                        )

        return sensitivity, specificity, accuracy, average_sensitivity_specificy

    def return_best_model(self, start_dir, database_signature):
        """
        Returns prediction model with the highest average cross-validation score

        RETURNS
            best_fit_model (obj.): object of class RandomForestClassifier with highest average cross-validation score
        """

        MODEL_FILENAME = self.training_taxon + str(database_signature) + '_assemblyQC.joblib.gz'
        COMPRESSION_TYPE = 'gzip'
        COMPRESSION_LEVEL = 3

        inclusion_data, exclusion_data = self.__curate_labels()
        print(inclusion_data.shape[0], exclusion_data.shape[0])
        list_models_and_scores = self.__generate_models(inclusion_data, exclusion_data)
        MODEL_INDEX = 0
        PERFORMANCE_METRICS_INDEX = 1
        AVERAGE_SCORE_INDEX = 3

        max_average_score = max(list_models_and_scores[i][PERFORMANCE_METRICS_INDEX][AVERAGE_SCORE_INDEX] for i in range(0, len(list_models_and_scores)))
        max_index = [list_models_and_scores[i][PERFORMANCE_METRICS_INDEX][AVERAGE_SCORE_INDEX] for i in range(0, len(list_models_and_scores))].index(max_average_score)
        print('Random Forest model {} has highest average sensitivity and specificity ({}). '.format(int(max_index+1), max_average_score), end='')
        print('Returning model {} as a python object.'.format(int(max_index+1)))
        best_fit_model = list_models_and_scores[max_index][MODEL_INDEX]

        joblib.dump(best_fit_model, os.path.join(start_dir, MODEL_FILENAME),
                    compress=(COMPRESSION_TYPE, COMPRESSION_LEVEL))

        model_object = joblib.load(os.path.join(start_dir, MODEL_FILENAME))

        return model_object
