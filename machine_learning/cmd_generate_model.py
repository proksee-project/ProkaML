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
import os
from pathlib import Path
from build_models import MachineLearningClassifier
import joblib

START_DIR = Path(__file__).resolve().parents[1]

filepath = '{}/well_represented_species_metadata_normalized.txt'.format(str(START_DIR))
dataframe = pd.read_table(filepath, low_memory=False)

classifier = MachineLearningClassifier(dataframe)
best_fit_model = classifier.return_best_model()

# Exporting best fitting model as a *.joblib python object
joblib.dump(best_fit_model,
            os.path.join(START_DIR, 'random_forest_n50_contigcount_l50_totlen_gccontent.joblib.gz'),
            compress=('gzip', 3))

dataframe_appended_probability = classifier.apply_model_to_database(best_fit_model)
dataframe_appended_probability.to_csv(os.path.join(START_DIR, 'well_represented_species_metadata_qc_probabilities.txt'),
                                      sep='\t', mode='w', index=False)
print('Assembly QC prediction probabilities appended to database')
