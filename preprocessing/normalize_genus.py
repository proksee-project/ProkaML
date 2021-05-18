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


class GenusNormalization():

    # Constants for normalized dataframe
    GENUS = 'Genus'
    SPECIES = 'Organism Name'

    # Constants for median attribute list indices
    N50 = 0
    NUM_CONTIGS = 1
    L50 = 2
    LENGTH = 3
    COVERAGE = 4
    GC_CONTENT = 5

    def __init__(self, species_median_dataframe, species_normalized_dataframe):

        self.species_median_dataframe = species_median_dataframe
        self.species_normalized_dataframe = species_normalized_dataframe
        self.median_normalization = MedianNormalization()
        self.genus_median_dict = {}

    def group_genus(self):

        TAXON = 'Species/Genus'
        GENUS_SPLIT_INDEX = 0
        SPLIT_CHAR = ' '

        self.species_median_dataframe[self.GENUS] = self.species_median_dataframe[TAXON].str.split(SPLIT_CHAR).str[GENUS_SPLIT_INDEX]
        genus_dataframe_group = self.species_median_dataframe.groupby(self.GENUS)
        genus_dataframe_list = [genus_dataframe_group.get_group(x) for x in genus_dataframe_group.groups]

        return genus_dataframe_list

    def get_genus_median_dictionary(self, dataframe, column_index):
        """
        Writes species specific median values of assembly attributes

        PARAMETERS
            dataframe (obj): object of class pandas.Dataframe with number of rows ranging from ten to several thousands,
            and 23 columns of different assembly attributes (str, int, float)

        POST
            species specific median attributes are written to MEDIAN_DATABASE_FILE
        """

        ROW_INDEX = 0
        genus = dataframe.iloc[ROW_INDEX, column_index]
        genus_median_logn50 = self.median_normalization.get_median_log_n50(dataframe)
        genus_median_logcontigcount = self.median_normalization.get_median_log_contigcount(dataframe)
        genus_median_logl50 = self.median_normalization.get_median_log_l50(dataframe)
        genus_median_logtotlen = self.median_normalization.get_median_log_totlen(dataframe)
        genus_median_logcoverage = self.median_normalization.get_median_log_coverage(dataframe)
        genus_median_gc_content = self.median_normalization.get_median_gc_content(dataframe)

        list_genus_median_attributes = [genus_median_logn50, genus_median_logcontigcount, genus_median_logl50,
                                        genus_median_logtotlen, genus_median_logcoverage, genus_median_gc_content]

        self.genus_median_dict[genus] = list_genus_median_attributes

        return list_genus_median_attributes

    def write_median_values(self, median_database_file):
        GENUS_COLUMN_INDEX = 7

        genus_dataframe_list = self.group_genus()
        for genus_dataframe in genus_dataframe_list:
            list_genus_median_attributes = self.get_genus_median_dictionary(genus_dataframe, GENUS_COLUMN_INDEX)
            median_database_file.write('{}\t'.format(genus_dataframe.iloc[0, GENUS_COLUMN_INDEX]))
            median_database_file.write('{}\t'.format(list_genus_median_attributes[self.L50]))
            median_database_file.write('{}\t'.format(list_genus_median_attributes[self.NUM_CONTIGS]))
            median_database_file.write('{}\t'.format(list_genus_median_attributes[self.L50]))
            median_database_file.write('{}\t'.format(list_genus_median_attributes[self.LENGTH]))
            median_database_file.write('{}\t'.format(list_genus_median_attributes[self.COVERAGE]))
            median_database_file.write('{}\n'.format(list_genus_median_attributes[self.GC_CONTENT]))


    def apply_genus_normalization(self):

        def get_genus_median_attribute(row, index):

            return self.genus_median_dict[row[self.GENUS]][index]

        logn50_normalized = 'logn50_normalized_genus'
        logcontigcount_normalized = 'logcontigcount_normalized_genus'
        logl50_normalized = 'logl50_normalized_genus'
        logtotlen_normalized = 'logtotlen_normalized_genus'
        logcoverage_normalized = 'logcoverage_normalized_genus'
        gccontent_normalized = 'gccontent_normalized_genus'

        self.species_normalized_dataframe[logn50_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.LOG_N50] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.N50), axis=1)

        self.species_normalized_dataframe[logcontigcount_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.LOG_NUM_CONTIGS] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.NUM_CONTIGS), axis=1)

        self.species_normalized_dataframe[logl50_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.LOG_L50] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.L50), axis=1)

        self.species_normalized_dataframe[logtotlen_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.LOG_LENGTH] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.LENGTH), axis=1)

        self.species_normalized_dataframe[logcoverage_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.LOG_COVERAGE] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.COVERAGE), axis=1)

        self.species_normalized_dataframe[gccontent_normalized] = \
            self.species_normalized_dataframe[self.median_normalization.GC_CONTENT_IMPUTED] - \
            self.species_normalized_dataframe.apply(lambda row: get_genus_median_attribute(row, self.GC_CONTENT), axis=1)

        return self.species_normalized_dataframe
