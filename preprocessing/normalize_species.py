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
from calculate_log_median import MedianNormalization


class SpeciesNormalization():
    """
    A class representing species specific normalization of assembly attributes

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
    """

    # Defining constants for column names
    SPECIES = 'Organism Name'
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
        self.median_normalization = MedianNormalization()

    def __group_species(self):
        """
        Subsets species with valid taxonomical groups

        RETURNS
            valid_species_dataframe_list (list): list of dataframes. Each dataframe has two-dimensional data structure, with
            number of rows ranging from ten to several thousands, and 16 columns of different assembly attributes (str, int, float)
        """

        species_dataframe_group = self.dataframe.groupby(self.SPECIES)
        species_dataframe_list = [species_dataframe_group.get_group(x) for x in species_dataframe_group.groups]

        return species_dataframe_list

    def __write_median_values(self, dataframe, median_database_file):
        """
        Writes species specific median values of assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        POST
            species specific median attributes are written to MEDIAN_DATABASE_FILE
        """

        ROW_INDEX = 0
        SPECIES_COLUMN_INDEX = 0

        median_database_file.write('{}\t'.format(dataframe.iloc[ROW_INDEX, SPECIES_COLUMN_INDEX]))
        median_database_file.write('{}\t'.format(self.median_normalization.get_median_log_n50(dataframe)))
        median_database_file.write('{}\t'.format(self.median_normalization.get_median_log_contigcount(dataframe)))
        median_database_file.write('{}\t'.format(self.median_normalization.get_median_log_l50(dataframe)))
        median_database_file.write('{}\t'.format(self.median_normalization.get_median_log_totlen(dataframe)))
        median_database_file.write('{}\t'.format(self.median_normalization.get_median_log_coverage(dataframe)))
        median_database_file.write('{}\n'.format(self.median_normalization.get_median_gc_content(dataframe)))

    def __normalize_dataframe(self, dataframe):
        """
        Calculates and appends species specific normalized assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        RETURNS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 29 columns of different assembly attributes (str, int, float)
        """

        logn50_normalized = 'logn50_normalized_species'
        logcontigcount_normalized = 'logcontigcount_normalized_species'
        logl50_normalized = 'logl50_normalized_species'
        logtotlen_normalized = 'logtotlen_normalized_species'
        logcoverage_normalized = 'logcoverage_normalized_species'
        gccontent_normalized = 'gccontent_normalized_species'

        dataframe[logn50_normalized] = dataframe[self.LOG_N50] - self.median_normalization.get_median_log_n50(dataframe)
        dataframe[logcontigcount_normalized] = dataframe[self.LOG_NUM_CONTIGS] - self.median_normalization.get_median_log_contigcount(dataframe)
        dataframe[logl50_normalized] = dataframe[self.LOG_L50] - self.median_normalization.get_median_log_l50(dataframe)
        dataframe[logtotlen_normalized] = dataframe[self.LOG_LENGTH] - self.median_normalization.get_median_log_totlen(dataframe)
        dataframe[logcoverage_normalized] = dataframe[self.LOG_COVERAGE] - self.median_normalization.get_median_log_coverage(dataframe)
        dataframe[gccontent_normalized] = dataframe[self.GC_CONTENT_IMPUTED] - self.median_normalization.get_median_gc_content(dataframe)

        return dataframe

    def execute(self, median_database_file):
        """
        Appends normalized attributes to entire database of assembly attributes

        POST
            Normalized assembly attributes for every assembly are written to NORMALIZED_DATABASE_FILE
        """

        species_aggregated_normalized_dataframe = []

        species_dataframe_list = self.__group_species()
        for species_dataframe in species_dataframe_list:
            logarithm_appended_dataframe = self.median_normalization.append_logarithm_imputations(species_dataframe)
            self.__write_median_values(logarithm_appended_dataframe, median_database_file)
            species_normalized_dataframe = self.__normalize_dataframe(logarithm_appended_dataframe)
            species_aggregated_normalized_dataframe.append(species_normalized_dataframe)

        species_aggregated_normalized_dataframe_concatenated = pd.concat(species_aggregated_normalized_dataframe)

        return species_aggregated_normalized_dataframe_concatenated
