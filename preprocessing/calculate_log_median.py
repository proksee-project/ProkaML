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


class MedianNormalization():
    """
    A class representing species specific normalization of assembly attributes

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
    """

    # Defining constants for column names
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


    def __init__(self):
        """
        Initializes the SpeciesNormalization class

        PARAMETERS
            dataframe (obj): dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
        """

        pass

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

        median_gc_content = dataframe[self.GC_CONTENT].median()

        return median_gc_content
