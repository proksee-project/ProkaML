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
import urllib.request
import urllib.error
import gzip
os.environ['OPENBLAS_NUM_THREADS'] = '1'
OUTPUT_DIRECTORY = 'additional_species_metadata'

class GCContentCalculate():
    """
    A class representing calculation of overall GC content of an assembly

    ATTRIBUTES
        email (str): User email address
        api_key (str): NCBI API key corresponding to user email address
        filename (str): A file containing metdata of assemblies (rows) and genomic attributes (columns) from NCBI
        output_dir (str): Directory for downloading contig assemblies
    """

    def __init__(self, email, api_key, dataframe):
        """
        Initializes the CalculateGCContent class

        PARAMETERS
            email (str): User email address
            api_key (str): NCBI API key corresponding to user email address
            filename (str): A file containing metdata of assemblies (rows) and genomic attributes (columns) from NCBI
            output_dir (str): Directory for downloading fasta contig assemblies
        """

        self.email = email
        self.api_key = api_key
        self.dataframe = dataframe

        if not os.path.exists(OUTPUT_DIRECTORY):
            os.mkdir(OUTPUT_DIRECTORY)

    def get_genbank_id_list(self):
        """
        Reads the metadata file

        RETURNS
            df (dataframe), genbank_id_list (list): tuple of dataframe and list of Genbank accession IDs
        """
        GENBANK_ID = 'Genbank Accession'
        genbank_id_list = self.dataframe[GENBANK_ID].to_list()

        return genbank_id_list

    def download_assembly(self, genbank_id):
        """
        Attempts to download genomic assembly using NCBI's API

        RETURNS
            assembly_file_path_local (str): local path to the downloaded assembly
        """

        # Constants for API queries
        DATABASE = 'assembly'
        ID_LIST = 'IdList'
        DOCUMENT_SUMMARY_SET = 'DocumentSummarySet'
        DOCUMENT_SUMMARY = 'DocumentSummary'
        INDEX = 0
        FTP_PATH_GENBANK = 'FtpPath_GenBank'
        VALIDATE = False
        ASSEMBLY_FILE_EXTENSION = '_genomic.fna.gz'

        try:
            # Program download link of assembly file
            handle = Entrez.esearch(db=DATABASE, term=genbank_id)
            record = Entrez.read(handle)
            handle2 = Entrez.esummary(db=DATABASE, id=record[ID_LIST])
            record2 = Entrez.read(handle2, validate=VALIDATE)
            genbank_dir = record2[DOCUMENT_SUMMARY_SET][DOCUMENT_SUMMARY][INDEX][FTP_PATH_GENBANK]
            assembly_file = os.path.basename(genbank_dir) + ASSEMBLY_FILE_EXTENSION
            assembly_file_link = os.path.join(genbank_dir, assembly_file)

            # Define path to download the assembly file
            assembly_file_path_local = os.path.join(OUTPUT_DIRECTORY, assembly_file)

            # Avoid duplicate downloading if assembly file exists
            if self.need_download(assembly_file_link, assembly_file_path_local):
                urllib.request.urlretrieve(assembly_file_link, assembly_file_path_local)

        except Exception:
            # Accounting for NCBI server issue or incorrect assembly file link
            assembly_file_path_local = 'NA'

        return assembly_file_path_local

    def need_download(self, url_link, assembly_file_path_local):
        """
        Checks for presence or absence of assembly file

        RETURNS
            bool : whether assembly file exists or not
        """

        url_request = urllib.request.Request(url_link, method='HEAD')
        url_open = urllib.request.urlopen(url_request, timeout=100)
        file_size = int(url_open.headers['Content-Length'])

        if (not os.path.exists(assembly_file_path_local) or
                file_size != os.path.getsize(assembly_file_path_local)):
            return True

        else:
            return False

    def calculate_gc_content(self, assembly_file_path_local):
        """
        Calculates GC content of an assembly as a fraction of the total assembly length

        RETURNS
            overall_gc_content (float): the overall GC content of an assembly
        """

        if assembly_file_path_local == 'NA':
            overall_gc_content = float('NaN')

        else:
            gc_content = 0
            full_length = 0
            try:
                with gzip.open(assembly_file_path_local, mode='rt') as open_file:
                    for line in open_file:
                        if not line.startswith('>'):
                            full_length += len(line.rstrip('\n'))
                            nucleotide_uppercase = line.rstrip('\n').upper()
                            for gc_pattern in 'GC':
                                gc_content += nucleotide_uppercase.count(gc_pattern)

                overall_gc_content = round(gc_content/full_length, 3)

            except Exception:
                # Accounting for zlib file opening errors
                overall_gc_content = float('NaN')

        return overall_gc_content

    def append_gc_content(self):
        """
        GC content of respective assemblies are appended as an additional column

        POST
            Creates file with additional column of GC content
        """

        Entrez.email = self.email
        Entrez.api_key = self.api_key

        genbank_id_list = self.get_genbank_id_list()

        # Iterate through list of Genbank IDs to download assembly file/s and calculate GC content
        gc_content_list = []
        for i in range(0, len(genbank_id_list)):
            assembly_file_path_local = self.download_assembly(genbank_id_list[i])
            print('{} assembly downloaded.'.format(genbank_id_list[i]), end='')

            # Calculate GC content and append to a list
            gc_content_list.append(self.calculate_gc_content(assembly_file_path_local))
            print('gc content calculated. ', end='')

            # Remove the downloaded file. This if block of code can be commented out/removed if 
            # the user wants to retain downloaded contig assemblies

            if assembly_file_path_local != 'NA':
                os.remove(assembly_file_path_local)
                print('{} assembly file removed'.format(genbank_id_list[i]))

        return gc_content_list
