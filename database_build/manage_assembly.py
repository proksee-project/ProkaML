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
from Bio import Entrez
import constants as const
import urllib.request
import urllib.error
import gzip
import re
os.environ['OPENBLAS_NUM_THREADS'] = '1'
SERVER_ERROR = 'NCBI server connection error'
MISMATCH_ERROR = 'esummary assembly link and actual link mismatch'
ZLIB_ERROR = 'zlib file opening error'


class AssemblyManager():
    """
    A class representing downloading and calculation of overall GC content of a genomic assembly

    ATTRIBUTES
        email (str): User email address
        api_key (str): NCBI API key corresponding to user email address
        dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
        species_log_file (str): Log file for every species cataloging assembly counts
    """

    def __init__(self, email, api_key, dataframe, species_log_file):
        """
        Initializes the AssemblyManager class

        PARAMETERS
            email (str): User email address
            api_key (str): NCBI API key corresponding to user email address
            dataframe (obj): species specific 2 dimensional dataframe with rows of assemblies and columns of genomic attributes
            species_log_file (str): Log file for every species cataloging assembly counts
        """

        self.email = email
        self.api_key = api_key
        self.dataframe = dataframe
        self.species_log_file = species_log_file

    def __get_genbank_id_list(self):
        """
        Reads the metadata file

        RETURNS
            genbank_id_list (list): list of Genbank accession IDs
        """

        genbank_id_list = self.dataframe[const.Metadata.METADATA_COLUMN_HEADERS[const.Metadata.METADATA_INDEX_GENBANK]].to_list()

        return genbank_id_list

    def download_assembly(self, genbank_id):
        """
        Attempts to download genomic assembly using NCBI's API

        RETURNS
            assembly_file_path_local (str): local path to the downloaded assembly
            error (str): error message while attempting to download assembly
        """

        ASSEMBLY_INDEX = 0
        error = const.FileFormat.EMPTY_STRING

        for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
            try:
                handle = Entrez.esearch(db=const.Assembly.ASSEMBLY_DATABASE, term=genbank_id)

            except Exception:
                assembly_file_path_local = const.FileFormat.NA

            else:
                record = Entrez.read(handle)
                handle2 = Entrez.esummary(db=const.Assembly.ASSEMBLY_DATABASE, id=record[const.Assembly.RECORD_IDLIST])
                record2 = Entrez.read(handle2, validate=const.Assembly.ESUMMARY_VALIDATE)
                genbank_dir = record2[const.Assembly.DOCUMENT_SUMMARY_SET][const.Assembly.DOCUMENT_SUMMARY][ASSEMBLY_INDEX][const.Assembly.FTP_PATH_GENBANK]
                if genbank_dir != const.FileFormat.EMPTY_STRING:
                    assembly_file = os.path.basename(genbank_dir) + const.Assembly.ASSEMBLY_FILE_EXTENSION
                    assembly_file_link = os.path.join(genbank_dir, assembly_file)
                    assembly_file_path_local = os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.ASSEMBLY_DOWNLOAD_DIR, assembly_file)
                else:
                    assembly_file_path_local = const.FileFormat.NA
                break

        if assembly_file_path_local == const.FileFormat.NA:
                error += SERVER_ERROR

        else:
            if self.__need_download(assembly_file_link, assembly_file_path_local):
                for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
                    try:
                        urllib.request.urlretrieve(assembly_file_link, assembly_file_path_local)

                    except Exception:
                        assembly_file_path_local = const.FileFormat.NA

                    else:
                        break

            if assembly_file_path_local == const.FileFormat.NA:
                error += MISMATCH_ERROR

        return assembly_file_path_local, error

    def __need_download(self, url_link, assembly_file_path_local):
        """
        Checks for presence or absence of assembly file

        RETURNS
            bool : whether assembly file exists or not
        """

        METHOD = 'HEAD'
        TIMEOUT = 100
        CONTENT_LENGTH = 'Content-Length'

        url_request = urllib.request.Request(url_link, method=METHOD)
        url_open = urllib.request.urlopen(url_request, timeout=TIMEOUT)
        file_size = int(url_open.headers[CONTENT_LENGTH])

        if (not os.path.exists(assembly_file_path_local) or
                file_size != os.path.getsize(assembly_file_path_local)):
            return True

        else:
            return False

    def __calculate_gc_content(self, assembly_file_path_local):
        """
        Calculates GC content of an assembly as a fraction of the total assembly length

        RETURNS
            overall_gc_content (float): the overall GC content of an assembly
            error (str): error message while unpacking assembly
        """

        ROUNDING_DIGITS = 3
        FASTA_FORMAT_START = '>'
        GC = 'GC'
        error = const.FileFormat.EMPTY_STRING
        gc_content = 0
        full_length = 0

        if assembly_file_path_local == const.FileFormat.NA:
            overall_gc_content = float('NaN')

        else:
            try:
                unzipped_assembly = gzip.open(assembly_file_path_local, mode='rt')

            except Exception:
                error += ZLIB_ERROR
                overall_gc_content = float('NaN')

            else:
                for line in unzipped_assembly:
                    if not line.startswith(FASTA_FORMAT_START):
                        full_length += len(line.rstrip('\n'))
                        nucleotide_uppercase = line.rstrip('\n').upper()
                        for gc_pattern in GC:
                            gc_content += nucleotide_uppercase.count(gc_pattern)

                overall_gc_content = round(gc_content/full_length, ROUNDING_DIGITS)

        return overall_gc_content, error

    def append_gc_content(self):
        """
        GC content of respective assemblies are appended as an additional column

        POST
            Creates file with additional column of GC content
        """

        Entrez.email = self.email
        Entrez.api_key = self.api_key
        genbank_id_list = self.__get_genbank_id_list()
        gc_content_list = []
        num_success = 0
        failed_id_server = []
        failed_id_mismatch = []
        failed_id_zlib = []

        for i in range(0, len(genbank_id_list)):
            assembly_file_path_local, download_error = self.download_assembly(genbank_id_list[i])

            if re.search(r'server', download_error):
                failed_id_server.append(genbank_id_list[i])
                print('{} assembly failed to download.'.format(genbank_id_list[i]))

            elif re.search(r'mismatch', download_error):
                failed_id_mismatch.append(genbank_id_list[i])
                print('{} assembly failed to download.'.format(genbank_id_list[i]))

            else:
                print('{} assembly downloaded, '.format(genbank_id_list[i]), end='')

            overall_gc_content, unzip_error = self.__calculate_gc_content(assembly_file_path_local)
            gc_content_list.append(overall_gc_content)

            if unzip_error:
                failed_id_zlib.append(genbank_id_list[i])
                print(' zlib error in reading assembly')
            elif not unzip_error and assembly_file_path_local != const.FileFormat.NA:
                num_success += 1
                print('gc content calculated.')

                # Remove the downloaded file. This if block of code can be commented out/removed if 
                # the user wants to retain downloaded contig assemblies
                os.remove(assembly_file_path_local)
                print('{} assembly file removed'.format(genbank_id_list[i]))

        return gc_content_list, num_success, failed_id_server, failed_id_mismatch, failed_id_zlib
