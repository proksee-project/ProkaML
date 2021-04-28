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


class OrganizeMetadata():
    """
    A class representing text uniformization of assembly methods and sequencing technologies in NCBI assembly database

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
    """

    def __init__(self, dataframe):
        """
        Initializes the OrganizeMetadata class

        PARAMETERS
            dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
        """

        self.dataframe = dataframe

    def organize_assembly_method(self):
        """
        Performs text cleaning of user submitted assembler methods by regular expressions based pattern
        identification and subsequent appropriate text replacements

        POST
            The 'Assembly Method' column data (str) is cleaned for clarity (str)
        """

        # Organizing assembly methods represented at least 10 times in the database
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Aa]5.*$'},
                                                {'Assembly Method': 'A5-miseq'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^A[Bb][Yy][Ss]{2}.*$'},
                                                {'Assembly Method': 'ABYSS'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Aa][Ll]{2}[Pp][Aa][Tt][Hh][Ss].*$'},
                                                {'Assembly Method': 'allpaths'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^AMO[Ss]cmp.*$'},
                                                {'Assembly Method': 'AMOScmp'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^ARGO.*$'},
                                                {'Assembly Method': 'ARGO'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Bb]owtie.*$'},
                                                {'Assembly Method': 'Bowtie'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^CA.*$'},
                                                {'Assembly Method': 'Celera'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Celera.*$'},
                                                {'Assembly Method': 'Celera'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Cc][Aa][Nn][Uu].*'},
                                                {'Assembly Method': 'Canu'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Cc][Ll][Cc].*$'},
                                                {'Assembly Method': 'CLC'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*DNA[Ss][Tt][Aa][Rr].*$'},
                                                {'Assembly Method': 'DNASTAR'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Seq[Mm]an.*$'},
                                                {'Assembly Method': 'DNASTAR'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Ee][Dd].*?[Nn][Aa].*$'},
                                                {'Assembly Method': 'Edena'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Ff][Aa][Ll][Cc][Oo][Nn].*$'},
                                                {'Assembly Method': 'Falcon'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Flye.*$'},
                                                {'Assembly Method': 'Flye'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Geneious.*$'},
                                                {'Assembly Method': 'Geneious'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Gg][Ss].*$'},
                                                {'Assembly Method': 'GS'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Roche\sGS.*$'},
                                                {'Assembly Method': 'GS'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Hh][Gg][Aa][Pp].*$'},
                                                {'Assembly Method': 'HGAP'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*Hierarchical\sGenome.*$'},
                                                {'Assembly Method': 'HGAP'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Ii][Dd][Bb][Aa].*$'},
                                                {'Assembly Method': 'IDBA_UD'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^INNU[Cc][Aa].*$'},
                                                {'Assembly Method': 'INNUca'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^MaSu[Rr]CA.*'},
                                                {'Assembly Method': 'MaSuRCA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Mm][Ee][Gg][Aa][Hh][Ii][Tt].*$'},
                                                {'Assembly Method': 'Megahit'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^MetaBAT.*$'},
                                                {'Assembly Method': 'MetaBAT'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Mm][Ii][Rr][Aa].*$'},
                                                {'Assembly Method': 'MIRA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Nn]ewbler.*$'},
                                                {'Assembly Method': 'Newbler'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^PATRIC.*$'},
                                                {'Assembly Method': 'PATRIC'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^pilon.*$'},
                                                {'Assembly Method': 'pilon'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Platanus.*$'},
                                                {'Assembly Method': 'Platanus'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^ProkaryoteAssembly.*$'},
                                                {'Assembly Method': 'ProkaryoteAssembly'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^QUAST.*$'},
                                                {'Assembly Method': 'QUAST'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Ray.*$'},
                                                {'Assembly Method': 'Ray'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Ss]hovill.*$'},
                                                {'Assembly Method': 'Shovill'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Ss][Kk][Ee][Ss][Aa].*$'},
                                                {'Assembly Method': 'SKESA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*SMRT.*$'},
                                                {'Assembly Method': 'SMRT'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^S[Oo][Aa][Pp].*$'},
                                                {'Assembly Method': 'SOAPdenovo'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Ss][Pp][Aa][Dd][Ee].*$'},
                                                {'Assembly Method': 'SPAdes'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Uu]ni[Cc]y.*$'},
                                                {'Assembly Method': 'unicycler'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*[Vv]el.+?t.*$'},
                                                {'Assembly Method': 'Velvet'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Ww]gs.*$'},
                                                {'Assembly Method': 'WgsAssembler'}, regex=True)

        # Replacing assembly methods with < 10 counts with NA
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Allora.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Arachne.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Artemis.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^As\s.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^BOI.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^breseq.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^BWA.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^CAP3.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^CISA.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Clover.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^custom.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^De.*?[Nn]ovo.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Galaxy.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^GeneStudio.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Illumina.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'.*in-house.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^IonTorrent.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Kbase.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Lasergene.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Mauve.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Microbe.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^[Mm]inia.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^MIX.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^MyPro.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^nanopolish.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^NextG[Ee][Nn]e.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^NovoAlign.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^other.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Parallel.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^PRJNA.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^prokka.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Shasta.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^SoftGenetics.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^SSPACE.*'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^TRIMMOMATIC.*'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Turing.*'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^Unknown.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({'Assembly Method': r'^VAAL.*$'},
                                                {'Assembly Method': 'NA'}, regex=True)

    def organize_sequencing_technology(self):
        """
        Performs text cleaning of user submitted sequencing technologies by regular expressions based pattern
        identification and subsequent appropriate text replacements

        POST
            The 'Sequencing Technology' column data (str) is cleaned for clarity (str)
        """

        # Textual organization of well known sequencing technologies
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'.*454.*$'},
                                                {'Sequencing Technology': '454'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^BGI.*$'},
                                                {'Sequencing Technology': 'BGI'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^Complete\sGenomics.*$'},
                                                {'Sequencing Technology': 'CompleteGenomics'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^DNB-Seq.*$'},
                                                {'Sequencing Technology': 'BGI'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[Ii].*?[Aa].*$'},
                                                {'Sequencing Technology': 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[HMhm][Ii][Ss][Ee][Qq].*$'},
                                                {'Sequencing Technology': 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[Nn]ext[Ss][Ee][Qq].*$'},
                                                {'Sequencing Technology': 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[Nn]ova[Ss][Ee][Qq].*$'},
                                                {'Sequencing Technology': 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^Nextera.*$'},
                                                {'Sequencing Technology': 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[Ii][Oo][Nn].*$'},
                                                {'Sequencing Technology': 'IonTorrent'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^Oxford\sNanopore.*$'},
                                                {'Sequencing Technology': 'OxfordNanopore'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'.*Min[Ii][Oo][Nn].*$'},
                                                {'Sequencing Technology': 'OxfordNanopore'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^[Pp]ac[Bb]io.*$'},
                                                {'Sequencing Technology': 'PacBio'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^Sanger.*$'},
                                                {'Sequencing Technology': 'Sanger'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'^Solexa.*$'},
                                                {'Sequencing Technology': 'Solexa'}, regex=True)
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'.*[Ss][Oo][Ll]i[Dd].*$'},
                                                {'Sequencing Technology': 'SOLiD'}, regex=True)

        # Replacing sequencing technologies that are not well known with 'NA'
        self.dataframe = self.dataframe.replace({'Sequencing Technology': r'.*Others.*$'},
                                                {'Sequencing Technology': 'NA'}, regex=True)

    def identify_long_reads(self):
        """
        Identifies assemblies associated with long reads

        RETURNS
            long_read_index (list) : list of integers that are row labels of the dataframe associated with long read
            assemblies
        """

        long_read_index = self.dataframe[(self.dataframe['Sequencing Technology'] == 'OxfordNanopore') |
                                         (self.dataframe['Sequencing Technology'] == 'PacBio') |
                                         (self.dataframe['Assembly Method'] == 'Canu') |
                                         (self.dataframe['Assembly Method'] == 'Falcon') |
                                         (self.dataframe['Assembly Method'] == 'Flye') |
                                         (self.dataframe['Assembly Method'] == 'HGAP') |
                                         (self.dataframe['Assembly Method'] == 'SMRT') |
                                         (self.dataframe['Assembly Method'] == 'pilon')].index

        return long_read_index

    def prenormalize_metadata(self):
        """
        Performs text cleaning by string replacement of assembly method and sequencing platform. 
        Excludes assembly rows associated with long reads.

        RETURNS
            self.dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
            ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
        """

        self.organize_assembly_method()
        self.organize_sequencing_technology()
        long_read_index = self.identify_long_reads()
        self.dataframe.drop(long_read_index, inplace=True)

        return self.dataframe
