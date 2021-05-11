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

import re

class CleanMetadata():
    """
    A class representing text uniformization of assembly methods and sequencing technologies in NCBI assembly database

    ATTRIBUTES
        dataframe (obj): An object of class pandas.Dataframe having a two-dimensional data structure with
        ~500,000 rows (first row contains headers) and 16 columns of different assembly attributes (str, int, float)
    """

    SPECIES = 'Organism Name'
    ASSEMBLY_METHOD = 'Assembly Method'
    SEQ_PLATFORM = 'Sequencing Technology'
    SPECIES_SPLIT_CHAR = ' '

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
        Re-writes similar assembly methods to a uniform assembly method name. For example: 'Spaded', 'SPAdes v 2.5',
        'SPAdes v 3.0', 'SPAdes assembler v 3.1.1', 'spades v 13.3' are all re-written to 'SPAdes'. A similar strategy
        is used for other assembly methods

        POST
            The 'Assembly Method' column (str) of metada is re-written so that similar assembly methods are same and uniform
        """

        # Organizing assembly methods represented at least 10 times in the database
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Aa]5.*$'},
                                                {self.ASSEMBLY_METHOD: 'A5-miseq'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^A[Bb][Yy][Ss]{2}.*$'},
                                                {self.ASSEMBLY_METHOD: 'ABYSS'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Aa][Ll]{2}[Pp][Aa][Tt][Hh][Ss].*$'},
                                                {self.ASSEMBLY_METHOD: 'allpaths'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^AMO[Ss]cmp.*$'},
                                                {self.ASSEMBLY_METHOD: 'AMOScmp'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^ARGO.*$'},
                                                {self.ASSEMBLY_METHOD: 'ARGO'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Bb]owtie.*$'},
                                                {self.ASSEMBLY_METHOD: 'Bowtie'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^CA.*$'},
                                                {self.ASSEMBLY_METHOD: 'Celera'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Celera.*$'},
                                                {self.ASSEMBLY_METHOD: 'Celera'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Cc][Aa][Nn][Uu].*'},
                                                {self.ASSEMBLY_METHOD: 'Canu'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Cc][Ll][Cc].*$'},
                                                {self.ASSEMBLY_METHOD: 'CLC'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*DNA[Ss][Tt][Aa][Rr].*$'},
                                                {self.ASSEMBLY_METHOD: 'DNASTAR'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Seq[Mm]an.*$'},
                                                {self.ASSEMBLY_METHOD: 'DNASTAR'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Ee][Dd].*?[Nn][Aa].*$'},
                                                {self.ASSEMBLY_METHOD: 'Edena'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Ff][Aa][Ll][Cc][Oo][Nn].*$'},
                                                {self.ASSEMBLY_METHOD: 'Falcon'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Flye.*$'},
                                                {self.ASSEMBLY_METHOD: 'Flye'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Geneious.*$'},
                                                {self.ASSEMBLY_METHOD: 'Geneious'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Gg][Ss].*$'},
                                                {self.ASSEMBLY_METHOD: 'GS'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Roche\sGS.*$'},
                                                {self.ASSEMBLY_METHOD: 'GS'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Hh][Gg][Aa][Pp].*$'},
                                                {self.ASSEMBLY_METHOD: 'HGAP'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*Hierarchical\sGenome.*$'},
                                                {self.ASSEMBLY_METHOD: 'HGAP'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Ii][Dd][Bb][Aa].*$'},
                                                {self.ASSEMBLY_METHOD: 'IDBA_UD'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^INNU[Cc][Aa].*$'},
                                                {self.ASSEMBLY_METHOD: 'INNUca'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^MaSu[Rr]CA.*'},
                                                {self.ASSEMBLY_METHOD: 'MaSuRCA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Mm][Ee][Gg][Aa][Hh][Ii][Tt].*$'},
                                                {self.ASSEMBLY_METHOD: 'Megahit'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^MetaBAT.*$'},
                                                {self.ASSEMBLY_METHOD: 'MetaBAT'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Mm][Ii][Rr][Aa].*$'},
                                                {self.ASSEMBLY_METHOD: 'MIRA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Nn]ewbler.*$'},
                                                {self.ASSEMBLY_METHOD: 'Newbler'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^PATRIC.*$'},
                                                {self.ASSEMBLY_METHOD: 'PATRIC'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^pilon.*$'},
                                                {self.ASSEMBLY_METHOD: 'pilon'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Platanus.*$'},
                                                {self.ASSEMBLY_METHOD: 'Platanus'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^ProkaryoteAssembly.*$'},
                                                {self.ASSEMBLY_METHOD: 'ProkaryoteAssembly'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^QUAST.*$'},
                                                {self.ASSEMBLY_METHOD: 'QUAST'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Ray.*$'},
                                                {self.ASSEMBLY_METHOD: 'Ray'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Ss]hovill.*$'},
                                                {self.ASSEMBLY_METHOD: 'Shovill'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Ss][Kk][Ee][Ss][Aa].*$'},
                                                {self.ASSEMBLY_METHOD: 'SKESA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*SMRT.*$'},
                                                {self.ASSEMBLY_METHOD: 'SMRT'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^S[Oo][Aa][Pp].*$'},
                                                {self.ASSEMBLY_METHOD: 'SOAPdenovo'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Ss][Pp][Aa][Dd][Ee].*$'},
                                                {self.ASSEMBLY_METHOD: 'SPAdes'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Uu]ni[Cc]y.*$'},
                                                {self.ASSEMBLY_METHOD: 'unicycler'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*[Vv]el.+?t.*$'},
                                                {self.ASSEMBLY_METHOD: 'Velvet'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Ww]gs.*$'},
                                                {self.ASSEMBLY_METHOD: 'WgsAssembler'}, regex=True)

        # Replacing assembly methods with < 10 counts with NA
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Allora.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Arachne.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Artemis.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^As\s.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^BOI.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^breseq.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^BWA.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^CAP3.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^CISA.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Clover.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^custom.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^De.*?[Nn]ovo.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Galaxy.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^GeneStudio.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Illumina.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'.*in-house.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^IonTorrent.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Kbase.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Lasergene.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Mauve.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Microbe.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^[Mm]inia.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^MIX.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^MyPro.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^nanopolish.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^NextG[Ee][Nn]e.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^NovoAlign.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^other.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Parallel.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^PRJNA.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^prokka.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Shasta.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^SoftGenetics.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^SSPACE.*'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^TRIMMOMATIC.*'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Turing.*'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^Unknown.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)
        self.dataframe = self.dataframe.replace({self.ASSEMBLY_METHOD: r'^VAAL.*$'},
                                                {self.ASSEMBLY_METHOD: 'NA'}, regex=True)

    def organize_sequencing_technology(self):
        """
        Re-writes similar sequencing technologies to a uniform sequencing technology name. For example: 'Illunina', 'Ilumina',
        'Illumina NovaSeq', 'Illumina NextSeq', 'ILLUMINA' are all re-written to 'Illumina'. A similar strategy
        is used for other sequencing technologies

        POST
            The 'Sequencing Technology' column (str) of metada is re-written so that similar sequencing technologies are same
            and uniform
        """

        # Textual organization of well known sequencing technologies
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'.*454.*$'},
                                                {self.SEQ_PLATFORM: '454'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^BGI.*$'},
                                                {self.SEQ_PLATFORM: 'BGI'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^Complete\sGenomics.*$'},
                                                {self.SEQ_PLATFORM: 'CompleteGenomics'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^DNB-Seq.*$'},
                                                {self.SEQ_PLATFORM: 'BGI'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[Ii].*?[Aa].*$'},
                                                {self.SEQ_PLATFORM: 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[HMhm][Ii][Ss][Ee][Qq].*$'},
                                                {self.SEQ_PLATFORM: 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[Nn]ext[Ss][Ee][Qq].*$'},
                                                {self.SEQ_PLATFORM: 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[Nn]ova[Ss][Ee][Qq].*$'},
                                                {self.SEQ_PLATFORM: 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^Nextera.*$'},
                                                {self.SEQ_PLATFORM: 'Illumina'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[Ii][Oo][Nn].*$'},
                                                {self.SEQ_PLATFORM: 'IonTorrent'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^Oxford\sNanopore.*$'},
                                                {self.SEQ_PLATFORM: 'OxfordNanopore'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'.*Min[Ii][Oo][Nn].*$'},
                                                {self.SEQ_PLATFORM: 'OxfordNanopore'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^[Pp]ac[Bb]io.*$'},
                                                {self.SEQ_PLATFORM: 'PacBio'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^Sanger.*$'},
                                                {self.SEQ_PLATFORM: 'Sanger'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'^Solexa.*$'},
                                                {self.SEQ_PLATFORM: 'Solexa'}, regex=True)
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'.*[Ss][Oo][Ll]i[Dd].*$'},
                                                {self.SEQ_PLATFORM: 'SOLiD'}, regex=True)

        # Replacing sequencing technologies that are not well known with 'NA'
        self.dataframe = self.dataframe.replace({self.SEQ_PLATFORM: r'.*Others.*$'},
                                                {self.SEQ_PLATFORM: 'NA'}, regex=True)

    def identify_long_reads(self):
        """
        Identifies assemblies associated with long reads

        RETURNS
            long_read_index (list) : list of integers that are row labels of the dataframe associated with long read
            assemblies
        """

        long_read_index = self.dataframe[(self.dataframe[self.SEQ_PLATFORM] == 'OxfordNanopore') |
                                         (self.dataframe[self.SEQ_PLATFORM] == 'PacBio') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'Canu') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'Falcon') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'Flye') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'HGAP') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'SMRT') |
                                         (self.dataframe[self.ASSEMBLY_METHOD] == 'pilon')].index

        return long_read_index

    def select_taxonomically_valid_species(self):
        """
        Subsets species with valid taxonomical groups

        RETURNS
            valid_species_dataframe_list (list): list of dataframes. Each dataframe has two-dimensional data structure, with
            number of rows ranging from ten to several thousands, and 16 columns of different assembly attributes (str, int, float)
        """

        REGEX_JOINING_CHAR = '|'
        SPECIES_EXPECTED_LENGTH = 2

        # List of regular expressions to exclude invalid species' taxonomy names, can be expanded
        EXCLUDE_PATTERNS = ['uncultured', 'metagenome']
        exclude_species_patterns = re.compile(REGEX_JOINING_CHAR.join(EXCLUDE_PATTERNS))

        # List of unexpected special characters in species' names
        UNEXPECTED_SPECIAL_CHARS = ['\[','\]']
        REPLACED_CHAR = ''
        unexpected_species_chars = re.compile(REGEX_JOINING_CHAR.join(UNEXPECTED_SPECIAL_CHARS))

        # Defining rules for filtering taxonomically valid species
        # Filtering condition 1: remove species with strings matching to EXCLUDE_PATTERNS
        filtering_condition1 = ~self.dataframe[self.SPECIES].str.contains(exclude_species_patterns)

        # Filtering condition 2: select only those species whose number of words equal SPECIES_EXPECTED_LENGTH
        filtering_condition2 = self.dataframe[self.SPECIES].str.split(self.SPECIES_SPLIT_CHAR).str.len() == SPECIES_EXPECTED_LENGTH
        self.dataframe = self.dataframe[filtering_condition1 & filtering_condition2]

        # Filtering condition 3: remove unexpected characters in species' names
        self.dataframe[self.SPECIES] = self.dataframe[self.SPECIES].str.replace(unexpected_species_chars, REPLACED_CHAR, regex=True)

    def assign_genus(self):
        GENUS_INDEX = 0
        self.dataframe = self.dataframe.assign(Genus=
                                               self.dataframe[self.SPECIES].str.split(self.SPECIES_SPLIT_CHAR).str[GENUS_INDEX])

        return self.dataframe

    def clean_metadata(self):
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
        self.select_taxonomically_valid_species()
        cleaned_dataframe = self.assign_genus()

        return cleaned_dataframe
