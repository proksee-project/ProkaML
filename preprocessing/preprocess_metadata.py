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

import numpy as np
import pandas as pd
import os
import re
from pathlib import Path

START_DIR = Path(__file__).resolve().parents[1]
MEDIAN_DATABASE_FILE = open(os.path.join(START_DIR, 'lineage_median_log_metrics.txt'), 'w')
MEDIAN_DATABASE_FILE.write('Species/Genus\tlogn50\tlogcontigcount\tlogl50\tlogtotlen\tlogcoverage\tgccontent\n')


class TaxonomicalNormalization():
    """
    A class representing species specific normalization of assembly attributes

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
    """

    # Defining constants for column names
    SPECIES = 'Organism Name'
    GENUS = 'Genus'
    N50 = 'ContigN50'
    LOG_N50 = 'logn50'
    NUM_CONTIGS = 'Contig count'
    LOG_NUM_CONTIGS = 'logcontigcount'
    L50 = 'ContigL50'
    LOG_L50 = 'logl50'
    LENGTH = 'Total length'
    LOG_LENGTH = 'logtotlen'
    COVERAGE = 'Genome Coverage'
    LOG_COVERAGE = 'logcoverage'
    GC_CONTENT = 'GCcontent'
    GC_CONTENT_IMPUTED = 'GCcontent_imputed'


    def __init__(self, dataframe):
        """
        Initializes the SpeciesNormalization class

        PARAMETERS
            dataframe (obj): dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
        """

        self.dataframe = dataframe

    def grouping_species(self):
        """
        Subsets species with valid taxonomical groups

        RETURNS
            valid_species_dataframe_list (list): list of dataframes. Each dataframe has two-dimensional data structure, with
            number of rows ranging from ten to several thousands, and 16 columns of different assembly attributes (str, int, float)
        """

        species_dataframe_group = self.dataframe.groupby(self.SPECIES)
        species_dataframe_list = [species_dataframe_group.get_group(x) for x in species_dataframe_group.groups]

        return species_dataframe_list

    def grouping_genus(self):

        genera_dataframe_group = self.dataframe.groupby(self.GENUS)
        genera_dataframe_list = [genera_dataframe_group.get_group(x) for x in genera_dataframe_group.groups]

        return genera_dataframe_list

    def append_logarithm_imputations(self, dataframe):
        """
        Appends logarithms (or imputations for missing data) of assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 16 columns of different assembly attributes (str, int, float)

        RETURNS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)
        """

        ROUNDING_DIGITS = 3
        ADJUSTED_COVERAGE = 'adj coverage'
        ZERO_COVERAGE_ADJUSTMENT = 0.1
        
        dataframe = dataframe.assign(logn50=round(np.log10(dataframe[self.N50]), ROUNDING_DIGITS))
        dataframe = dataframe.assign(logcontigcount=round(np.log10(dataframe[self.NUM_CONTIGS]), ROUNDING_DIGITS))
        dataframe = dataframe.assign(logl50=round(np.log10(dataframe[self.L50]), ROUNDING_DIGITS))
        dataframe = dataframe.assign(logtotlen=round(np.log10(dataframe[self.LENGTH]), ROUNDING_DIGITS))

        # Adding column ADJUSTED_COVERAGE to allow logarithm calculations of assemblies with zero coverage
        dataframe.loc[dataframe[self.COVERAGE] > 0, ADJUSTED_COVERAGE] = dataframe[self.COVERAGE]
        dataframe.loc[dataframe[self.COVERAGE] == 0, ADJUSTED_COVERAGE] = ZERO_COVERAGE_ADJUSTMENT
        dataframe = dataframe.assign(logcoverage=round(np.log10(dataframe[ADJUSTED_COVERAGE]), ROUNDING_DIGITS))

        # Add column GC_CONTENT_IMPUTED to replace missing GC values with species specific medians
        dataframe.loc[dataframe[self.GC_CONTENT].notnull(), self.GC_CONTENT_IMPUTED] = dataframe[self.GC_CONTENT]
        dataframe.loc[dataframe[self.GC_CONTENT].isnull(), self.GC_CONTENT_IMPUTED] = dataframe[self.GC_CONTENT].median()

        return dataframe

    def get_median_log_n50(self, dataframe):
        """
        Calculates median of log(n50)

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_n50 (float): species specific median of log(n50)
        """

        median_log_n50 = dataframe[self.LOG_N50].median()

        return median_log_n50

    def get_median_log_contigcount(self, dataframe):
        """
        Calculates median of log(number of contigs)

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_contigcount (float): species specific median of log(number of contigs)
        """

        median_log_contigcount = dataframe[self.LOG_NUM_CONTIGS].median()

        return median_log_contigcount

    def get_median_log_l50(self, dataframe):
        """
        Calculates median of log(l50)

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_l50 (float): species specific median of log(l50)
        """

        median_log_l50 = dataframe[self.LOG_L50].median()

        return median_log_l50

    def get_median_log_totlen(self, dataframe):
        """
        Calculates median of log(assembly length)

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_totlen (float): species specific median of log(assembly length)
        """

        median_log_totlen = dataframe[self.LOG_LENGTH].median()

        return median_log_totlen

    def get_median_log_coverage(self, dataframe):
        """
        Calculates median of log(genome coverage)

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_coverage (float): species specific median of log(genome coverage)
        """

        median_log_coverage = dataframe[self.LOG_COVERAGE].median()

        return median_log_coverage

    def get_median_gc_content(self, dataframe):
        """
        Calculates median of overall gc content

        ATTRIBUTES
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            median_log_coverage (float): species specific median of gc content
        """

        median_gc_content = dataframe[self.GC_CONTENT_IMPUTED].median()

        return median_gc_content

    def write_median_values(self, dataframe, column_index):
        """
        Writes species specific median values of assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        POST
            species specific median attributes are written to MEDIAN_DATABASE_FILE
        """

        MEDIAN_DATABASE_FILE.write('{}\t'.format(dataframe.iloc[0, column_index]))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_n50(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_contigcount(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_l50(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_totlen(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_log_coverage(dataframe)))
        MEDIAN_DATABASE_FILE.write('{}\t'.format(self.get_median_gc_content(dataframe)))
        MEDIAN_DATABASE_FILE.write('\n')

    def normalize_dataframe(self, dataframe, taxonomy_resolution):
        """
        Calculates and appends species specific normalized assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 29 columns of different assembly attributes (str, int, float)
        """

        logn50_normalized = 'logn50_normalized_' + str(taxonomy_resolution)
        logcontigcount_normalized = 'logcontigcount_normalized_' + str(taxonomy_resolution)
        logl50_normalized = 'logl50_normalized_' + str(taxonomy_resolution)
        logtotlen_normalized = 'logtotlen_normalized_' + str(taxonomy_resolution)
        logcoverage_normalized = 'logcoverage_normalized_' + str(taxonomy_resolution)
        gccontent_normalized = 'gccontent_normalized_' + str(taxonomy_resolution)

        dataframe[logn50_normalized] = dataframe[self.LOG_N50] - self.get_median_log_n50(dataframe)
        dataframe[logcontigcount_normalized] = dataframe[self.LOG_NUM_CONTIGS] - self.get_median_log_contigcount(dataframe)
        dataframe[logl50_normalized] = dataframe[self.LOG_L50] - self.get_median_log_l50(dataframe)
        dataframe[logtotlen_normalized] = dataframe[self.LOG_LENGTH] - self.get_median_log_totlen(dataframe)
        dataframe[logcoverage_normalized] = dataframe[self.LOG_COVERAGE] - self.get_median_log_coverage(dataframe)
        dataframe[gccontent_normalized] = dataframe[self.GC_CONTENT_IMPUTED] - self.get_median_gc_content(dataframe)

        return dataframe

    def apply_species_normalization(self):
        """
        Appends normalized attributes to entire database of assembly attributes

        POST
            Normalized assembly attributes for every assembly are written to NORMALIZED_DATABASE_FILE
        """

        SPECIES_TAXONOMY = 'species'
        SPECIES_COLUMN_INDEX = 0
        species_aggregated_normalized_dataframe = []

        species_dataframe_list = self.grouping_species()
        for species_dataframe in species_dataframe_list:
            logarithm_appended_dataframe = self.append_logarithm_imputations(species_dataframe)
            self.write_median_values(logarithm_appended_dataframe, SPECIES_COLUMN_INDEX)
            species_normalized_dataframe = self.normalize_dataframe(logarithm_appended_dataframe, SPECIES_TAXONOMY)
            species_aggregated_normalized_dataframe.append(species_normalized_dataframe)

        species_aggregated_normalized_dataframe_concatenated = pd.concat(species_aggregated_normalized_dataframe)

        return species_aggregated_normalized_dataframe_concatenated

    def apply_genus_normalization(self):

        GENUS_TAXONOMY = 'genus'
        GENUS_COLUMN_INDEX = 16
        genera_normalized_dataframe = []

        genera_dataframe_list = self.grouping_genus()
        for genus_dataframe in genera_dataframe_list:
            logarithm_appended_dataframe = self.append_logarithm_imputations(genus_dataframe)
            self.write_median_values(logarithm_appended_dataframe, GENUS_COLUMN_INDEX)
            genus_normalized_dataframe = self.normalize_dataframe(logarithm_appended_dataframe, GENUS_TAXONOMY)
            genera_normalized_dataframe.append(genus_normalized_dataframe)

        genera_aggregated_normalized_dataframe_concatenated = pd.concat(genera_normalized_dataframe)

        return genera_aggregated_normalized_dataframe_concatenated
