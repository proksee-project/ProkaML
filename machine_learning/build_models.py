'''
Copyright:
University of Manitoba & National Microbiology Laboratory, Canada, 2020
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

# Defining instance of RandomForestClassifier with parameters
rfc = RandomForestClassifier(n_estimators=100,
                             criterion="gini",
                             max_depth=None,
                             max_features="auto",
                             random_state=42)


class MachineLearningClassifier():
    """
    A class for generating and evaluating machine learning models from assembly QC data

    ATTRIBUTES
        dataframe (obj): the dataframe object of assembly metadata
    """

    def __init__(self, dataframe):
        """
        Initializes the MachineLearningClassifier class

        PARAMETERS
            dataframe (obj): the dataframe object of assembly metadata
        """

        self.dataframe = dataframe

    def curate_labels(self):
        """
        Curates labels based on assembly inclusion or exclusion in NCBI RefSeq database

        RETURNS
           inclusion_data, exclusion_data : tuple of dataframes (obj.) for inclusion and exclusion datasets
        """

        inclusion_label = 1
        exclusion_label = 2

        RefSeq_included_data = self.dataframe[self.dataframe['Refseq Accession'].notnull()]
        RefSeq_included_data = RefSeq_included_data.assign(label=inclusion_label)

        RefSeq_excluded_data = self.dataframe[self.dataframe['Refseq Accession'].isnull()]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_data[~RefSeq_excluded_data['Refseq Exclusion Reason'].
                                                             str.contains(r'multi-isolate|\[\]', regex=True)]
        RefSeq_excluded_multi_isolate = RefSeq_excluded_multi_isolate.assign(label=exclusion_label)

        # Creating dataframes of inclusion and exclusion datasets with columns of interest
        inclusion_data = RefSeq_included_data[['logn50_normalized', 'logcontigcount_normalized', 'logl50_normalized',
                                               'logtotlen_normalized', 'gccontent_normalized', 'label']]
        exclusion_data = RefSeq_excluded_multi_isolate[['logn50_normalized', 'logcontigcount_normalized',
                                                        'logl50_normalized', 'logtotlen_normalized',
                                                        'gccontent_normalized', 'label']]

        return inclusion_data, exclusion_data

    def generate_models(self, inclusion_data, exclusion_data, num_iterations):
        """
        Generates random forests models and calculates cross-validation scores

        PARAMETERS
            inclusion_data (obj.): dataframe of NCBI RefSeq included assemblies
            exclusion_data (obj.): dataframe of NCBI RefSeq excluded assemblies
            num_iterations (int): number of random forests models to be generated

        RETURNS
            list_models_and_scores (list): list of prediction models with averaged cross-validation scores
        """

        list_models_and_scores = []
        for i in range(0, num_iterations):

            # Randomized balanced subsamping of inclusion data
            inclusion_data_balanced = inclusion_data.sample(n=exclusion_data.shape[0])
            integrated_data = pd.concat([inclusion_data_balanced, exclusion_data], ignore_index=True)

            # Defining feature space and label vector
            feature_space = integrated_data.drop('label', axis=1)
            label_vector = integrated_data['label']

            # Fitting random forests model to the training data
            rfc.fit(feature_space, label_vector)
            print('Fitting random forest model {}. Evaluating cross-validation score'.format(int(i+1)))

            # Evaluating models by 10 fold cross-validation
            n_fold_cv = 10
            model_average_cv_score = np.average(cross_val_score(rfc, feature_space, label_vector, cv=n_fold_cv))

            list_models_and_scores.append([rfc, model_average_cv_score])

        return list_models_and_scores

    def return_best_model(self):
        """
        Returns prediction model with the highest average cross-validation score

        RETURNS
            best_fit_model (obj.): random forests model object with highest average cross-validation score
        """

        inclusion_data, exclusion_data = self.curate_labels()
        num_iterations = 10
        list_models_and_scores = self.generate_models(inclusion_data, exclusion_data, num_iterations)

        max_score = max(list_models_and_scores[i][1] for i in range(0, len(list_models_and_scores)))
        max_index = [list_models_and_scores[i][1] for i in range(0, len(list_models_and_scores))].index(max_score)
        print('Random Forest model {} has highest average cross-validation score. '.format(int(max_index+1)), end='')
        print('Returning model {} as a python object.'.format(int(max_index+1)))
        best_fit_model = list_models_and_scores[max_index][0]

        return best_fit_model

    def apply_model_to_database(self, best_fit_model):
        """
        Applies best fitting model to entire database

        PARAMETERS
            best_fit_model (obj.): random forests model object with highest average cross-validation score

        RETURNS
            self.dataframe (obj.): database object with added column of prediction probability
        """

        database_dataframe = self.dataframe[['logn50_normalized', 'logcontigcount_normalized', 'logl50_normalized',
                                             'logtotlen_normalized', 'gccontent_normalized']]

        # Predicting assembly QC probability and assigning prediction probability column to database
        self.dataframe['RefSeq_predict_prob'] = best_fit_model.predict_proba(database_dataframe)[:, 0]

        return self.dataframe
