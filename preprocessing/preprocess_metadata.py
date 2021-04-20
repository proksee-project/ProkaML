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
import os
import sys
import re 
from pathlib import Path

START_DIR = Path(__file__).resolve().parents[1]
MEDIAN_DATABASE_FILE = open(os.path.join(START_DIR, 'species_median_log_metrics.txt'), 'w')
MEDIAN_DATABASE_FILE.write('Species\tlogn50\tlogcontigcount\tlogl50\tlogtotlen\tlogcoverage\tgccontent\n')
NORMALIZED_DATABASE_FILE = open(os.path.join(START_DIR, 'well_represented_species_metadata_normalized.txt'), 'w')

class SpeciesNormalization():
    """
    A class representing species specific normalization of assembly attributes

    ATTRIBUTES
        dataframe (obj): the dataframe object of assembly metadata
    """

    def __init__(self, dataframe):
        """
        Initializes the SpeciesNormalization class

        PARAMETERS
            dataframe (obj): the dataframe object of assembly metadata
        """

        self.dataframe = dataframe

    def subset_valid_species(self):
        """
        Subsets species with valid taxonomical groups

        RETURNS
            valid_species_dataframe_list (list): list of dataframes of (taxonomically valid) species
            containing assembly attributes
        """

        # List of regular expressions to exclude invalid species' taxonomy names, can be expanded
        exclude_species_patterns = ['uncultured', 
                                    'metagenome',
                                    '^bacterium$'
                                    ]
        exclude_species = re.compile("|".join(exclude_species_patterns))
        valid_species_dataframe = self.dataframe[~self.dataframe['Organism Name'].str.contains(exclude_species)]
        valid_species_group = valid_species_dataframe.groupby('Organism Name')
        valid_species_dataframe_list = [valid_species_group.get_group(x) for x in valid_species_group.groups]

        return valid_species_dataframe_list

    def append_logarithm_imputations(self, dataframe):
        """
        Appends logarithms (or imputations for missing data) of assembly attributes

        PARAMETERS
            dataframe (obj): dataframe object of species specific attributes

        RETURNS
            dataframe (obj): dataframe with appended columns of species specific logarithm calculated attributes
        """

        dataframe = dataframe.assign(logn50=round(np.log10(dataframe['ContigN50']),3))
        dataframe = dataframe.assign(logcontigcount=round(np.log10(dataframe['Contig count']),3))
        dataframe = dataframe.assign(logl50=round(np.log10(dataframe['ContigL50']),3))
        dataframe = dataframe.assign(logtotlen=round(np.log10(dataframe['Total length']),3))

        # Adding column 'adj coverage' to allow logarithm calculations of assemblies with zero coverage
        dataframe.loc[dataframe['Genome Coverage'] > 0, 'adj coverage'] = dataframe['Genome Coverage']
        dataframe.loc[dataframe['Genome Coverage'] == 0, 'adj coverage'] = 0.1
        dataframe = dataframe.assign(logcoverage=round(np.log10(dataframe['adj coverage']),3))

        # Add column 'GCcontent_imputed' to replace missing GC values with species specific medians
        dataframe.loc[dataframe['GCcontent'].notnull(), 'GCcontent_imputed'] = dataframe['GCcontent']
        dataframe.loc[dataframe['GCcontent'].isnull(), 'GCcontent_imputed'] = dataframe['GCcontent'].median()

        return dataframe

    def get_median_log_n50(self, dataframe):
        """
        Calculates median of log(n50)

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_n50 (float): species specific median of log(n50)
        """

        median_log_n50 = dataframe['logn50'].median()

        return median_log_n50

    def get_median_log_contigcount(self, dataframe):
        """
        Calculates median of log(number of contigs)

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_contigcount (float): species specific median of log(number of contigs)
        """

        median_log_contigcount = dataframe['logcontigcount'].median()

        return median_log_contigcount

    def get_median_log_l50(self, dataframe):
        """
        Calculates median of log(l50)

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_l50 (float): species specific median of log(l50)
        """

        median_log_l50 = dataframe['logl50'].median()

        return median_log_l50

    def get_median_log_totlen(self, dataframe):
        """
        Calculates median of log(assembly length)

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_totlen (float): species specific median of log(assembly length)
        """

        median_log_totlen = dataframe['logtotlen'].median()

        return median_log_totlen

    def get_median_log_coverage(self, dataframe):
        """
        Calculates median of log(genome coverage)

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_coverage (float): species specific median of log(genome coverage)
        """

        median_log_coverage = dataframe['logcoverage'].median()

        return median_log_coverage

    def get_median_gc_content(self, dataframe):
        """
        Calculates median of overall gc content

        ATTRIBUTES
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            median_log_coverage (float): species specific median of gc content
        """

        median_gc_content = dataframe['GCcontent_imputed'].median()

        return median_gc_content

    def write_median_values(self, dataframe):
        """
        Writes species specific median values of assembly attributes

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        POST
            species specific median attributes are written to MEDIAN_DATABASE_FILE
        """

        MEDIAN_DATABASE_FILE.write('{}\t'.format(dataframe.iloc[0, 0]))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_n50(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_contigcount(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_l50(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_totlen(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_coverage(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_gc_content(dataframe)))
        MEDIAN_DATABASE_FILE.write('\n')

    def normalize_dataframe(self, dataframe):
        """
        Calculates and appends species specific normalized assembly attributes

        PARAMETERS
            dataframe (obj): A species specific dataframe of assembly attributes

        RETURNS
            dataframe (obj): dataframe with appended columns of normalized attributes
        """

        dataframe = dataframe.assign(logn50_normalized = dataframe['logn50']- self.get_median_log_n50(dataframe))
        dataframe = dataframe.assign(logcontigcount_normalized = dataframe['logcontigcount']- self.get_median_log_contigcount(dataframe))
        dataframe = dataframe.assign(logl50_normalized = dataframe['logl50']- self.get_median_log_l50(dataframe))
        dataframe = dataframe.assign(logtotlen_normalized = dataframe['logtotlen']- self.get_median_log_totlen(dataframe))
        dataframe = dataframe.assign(logcoverage_normalized = dataframe['logcoverage']- self.get_median_log_coverage(dataframe))
        dataframe = dataframe.assign(gccontent_normalized = dataframe['GCcontent_imputed']- self.get_median_gc_content(dataframe))

        return dataframe

    def apply_normalization_to_database(self):
        """
        Appends normalized attributes to entire database of assembly attributes

        POST
            Normalized assembly attributes for every assembly are written to NORMALIZED_DATABASE_FILE
        """

        valid_species_dataframe_list = self.subset_valid_species()
        species_count = 0

        for individual_species_dataframe in valid_species_dataframe_list:
            appended_species_dataframe = self.append_logarithm_imputations(individual_species_dataframe)
            self.write_median_values(appended_species_dataframe)
            species_normalized_dataframe = self.normalize_dataframe(appended_species_dataframe)

            # Column headers are written for first instance only
            if species_count == 0:
                header = True
            else:
                header = False

            species_normalized_dataframe.to_csv(NORMALIZED_DATABASE_FILE, sep='\t', mode='a', header=header, index=False)
            species_count += 1
