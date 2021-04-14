import pandas as pd
import numpy as np
import re
import os
from pathlib import Path 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score
from sklearn.impute import SimpleImputer
import joblib

rfc = RandomForestClassifier()

START_DIR = Path(__file__).resolve().parents[1]
filepath = '{}/preprocessing/well_represented_species_metadata_normalized.txt'.format(str(START_DIR))

#Creating refseq inclusion and exclusion (minus multi-isolates) dataframes
df = pd.read_table(filepath)
refseqY_df = df[df['Refseq Accession'].notnull()]
refseqN_df = df[df['Refseq Accession'].isnull()]
refseqY_df = refseqY_df.assign(datalabel=1)
refseqN_further = refseqN_df[~refseqN_df['Refseq Exclusion Reason'].str.contains('multi-isolate|\[\]')]
refseqN_further = refseqN_further.assign(datalabel=2)

#Creating dataframes of inclusion and exclusion sets with columns of interest
inclusion_data = refseqY_df[['logn50_normalized','logcontigcount_normalized','logl50_normalized',
                            'logtotlen_normalized','gccontent_normalized','datalabel']]
exclusion_data = refseqN_further[['logn50_normalized','logcontigcount_normalized','logl50_normalized',
                                 'logtotlen_normalized','gccontent_normalized','datalabel']]

integrated_frame = pd.concat([inclusion_data, exclusion_data], ignore_index=True)

def iterating_random_forest_models():
    model_and_score = []
    for i in range(0,10):
        #Subsamping inclusion frame to equalize exclusion frame
        inclusion_data_balanced = inclusion_data.sample(n=exclusion_data.shape[0])
        integrated_data = pd.concat([inclusion_data_balanced, exclusion_data], ignore_index=True)

        #Defining data matrix and label vector
        X_train = integrated_data.drop('datalabel', axis=1)
        Y_train = integrated_data['datalabel']

        #Some GCnormalized values are missing, so we perform median imputation
        imputed = SimpleImputer(missing_values=np.nan, strategy='median')
        imputed = imputed.fit(X_train)
        X_train_imputed = imputed.transform(X_train)

        #training and importing random forest model
        rfc.fit(X_train_imputed, Y_train)
        pred_prob = rfc.predict_proba(X_train_imputed)
        auc_score = roc_auc_score(Y_train, pred_prob[:,1])
        model_and_score.append([rfc, auc_score])
    
    return model_and_score

def evaluate_best_model(model_and_score):
    auc_list = []
    for i in range(0, len(model_and_score)):
        auc_list.append(model_and_score[i][1])
    max_auc = max(auc_list)
    max_index = auc_list.index(max_auc)

    return model_and_score[max_index][0]

iterated_models_and_auc = iterating_random_forest_models()
best_fit_model = evaluate_best_model(iterated_models_and_auc)
joblib.dump(best_fit_model, "random_forest_n50_contigcount_l50_totlen_gccontent.joblib")

X_test = integrated_frame.drop('datalabel', axis=1)
print(best_fit_model.predict_proba(X_test))