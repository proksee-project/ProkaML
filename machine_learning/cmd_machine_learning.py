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
import re
import os
from pathlib import Path 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.impute import SimpleImputer
import joblib

START_DIR = Path(__file__).resolve().parents[1]

def label_curation(filepath):

    #Creating refseq inclusion and exclusion (minus multi-isolates) dataframes
    df = pd.read_table(filepath)
    refseqY_df = df[df['Refseq Accession'].notnull()]
    refseqY_df = refseqY_df.assign(datalabel=1)

    refseqN_df = df[df['Refseq Accession'].isnull()]
    refseqN_further = refseqN_df[~refseqN_df['Refseq Exclusion Reason'].str.contains('multi-isolate|\[\]')]
    refseqN_further = refseqN_further.assign(datalabel=2)

    #Creating dataframes of inclusion and exclusion sets with columns of interest
    inclusion_data = refseqY_df[['logn50_normalized','logcontigcount_normalized','logl50_normalized',
                                'logtotlen_normalized','gccontent_normalized','datalabel']]
    exclusion_data = refseqN_further[['logn50_normalized','logcontigcount_normalized','logl50_normalized',
                                    'logtotlen_normalized','gccontent_normalized','datalabel']]

    return inclusion_data, exclusion_data

def iterating_random_forest_models(inclusion_data, exclusion_data):

    classifier_and_score = []
    number_of_iterations = 10
    for i in range(0,number_of_iterations):

        #Ramdom subsamping of inclusion data for equalizing with exclusion data
        inclusion_data_balanced = inclusion_data.sample(n=exclusion_data.shape[0])
        integrated_data = pd.concat([inclusion_data_balanced, exclusion_data], ignore_index=True)

        #Defining feature matrix and labels
        X_train = integrated_data.drop('datalabel', axis=1)
        Y_train = integrated_data['datalabel']

        #Some GCnormalized values are missing, so we perform median imputation
        imputed = SimpleImputer(missing_values=np.nan, strategy='median')
        imputed = imputed.fit(X_train)
        X_train_imputed = imputed.transform(X_train)

        #Evaluating models by 10 fold cross-validation
        model_average_score = np.average(cross_val_score(RandomForestClassifier(), X_train_imputed, Y_train, cv=10))

        # Fitting model to the training data
        RandomForestClassifier().fit(X_train_imputed, Y_train)

        classifier_and_score.append([RandomForestClassifier(), model_average_score])
    
    return classifier_and_score

def evaluate_best_model(classifier_and_score):

    # Returning classifier model with highest average score
    average_score_list = []
    for i in range(0, len(classifier_and_score)):
        average_score_list.append(classifier_and_score[i][1])
    max_score = max(average_score_list)
    max_index = average_score_list.index(max_score)

    return classifier_and_score[max_index][0]

def main():

    filepath = '{}/well_represented_species_metadata_normalized.txt'.format(str(START_DIR))
    inclusion, exclusion = label_curation(filepath)
    classifiers_and_scores = iterating_random_forest_models(inclusion, exclusion)
    best_fit_model = evaluate_best_model(classifiers_and_scores)
    joblib.dump(best_fit_model, "random_forest_n50_contigcount_l50_totlen_gccontent.joblib")

if __name__ == "__main__":
    main()
