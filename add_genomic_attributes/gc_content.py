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

import os
import pandas as pd
from Bio import Entrez
import urllib.request
import urllib.error
import gzip
os.environ['OPENBLAS_NUM_THREADS'] = '1'


class CalculateGCContent():
    """
    A class representing calculation of overall GC content of an assembly

    ATTRIBUTES
        email (str): User email address
        api_key (str): NCBI API key corresponding to user email address
        filename (str): A file containing metdata of assemblies (rows) and genomic attributes (columns) from NCBI
        output_dir (str): Directory for downloading contig assemblies
    """

    def __init__(self, email, api_key, filename, output_dir):
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
        self.filename = filename
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    def read_file(self):
        """
        Reads the metadata file

        RETURNS
            df (dataframe), genbank_id_list (list): tuple of dataframe and list of Genbank accession IDs
        """

        df = pd.read_csv(self.filename, sep='\t', keep_default_na=False)
        genbank_id_list = df['Genbank Accession'].to_list()

        return df, genbank_id_list

    def download_assembly(self, genbank_id):
        """
        Programs the download link of fasta assembly and attempts to download it

        RETURNS
            assembly_file_path_local (str): local path to the downloaded fasta assembly
        """

        try:
            # Program download link of assembly file
            handle = Entrez.esearch(db="assembly", term=genbank_id)
            record = Entrez.read(handle)
            handle2 = Entrez.esummary(db="assembly", id=record['IdList'])
            record2 = Entrez.read(handle2, validate=False)
            genbank_dir = record2['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            assembly_file = os.path.basename(genbank_dir) + '_genomic.fna.gz'
            assembly_file_link = os.path.join(genbank_dir, assembly_file)

            # Define path to download the assembly file
            assembly_file_path_local = os.path.join(self.output_dir, assembly_file)

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
            overall_gc_content = 'NA'

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
                overall_gc_content = 'NA'

        return overall_gc_content

    def append_gc_content(self):
        Entrez.email = self.email
        Entrez.api_key = self.api_key

        dataframe, genbank_id_list = self.read_file()

        # Iterate through list of Genbank IDs to download assembly file/s and calculate GC content
        gc_content_list = []
        for i in range(0, len(genbank_id_list)):
            assembly_file_path_local = self.download_assembly(genbank_id_list[i])
            print('{} assembly downloaded.'.format(genbank_id_list[i]), end='')

            # Calculate GC content and append to a list
            gc_content_list.append(self.calculate_gc_content(assembly_file_path_local))
            print('gc content calculated.', sep='')

            # Remove the downloaded file. This try except can be commented out/removed if assemblies do not
            # consume much hard drive space
            try:
                os.remove(assembly_file_path_local)
                print('{} assembly file removed'.format(genbank_id_list[i]))

            except FileNotFoundError:
                pass

        # Write the list of GC content as a separate column to output file
        dataframe['GCcontent'] = gc_content_list
        output_file = self.filename.split('.')[0] + '_added_attributes.txt'
        dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
        dataframe.to_csv(output_file, sep='\t', mode='w', index=False)
