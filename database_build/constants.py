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

from Bio import Entrez
import os

Entrez.max_tries = 1
Entrez.sleep_between_tries = 1
API_QUERY_ATTEMPT_START = 1
API_QUERY_ATTEMPT_END = 4


class FileFormat():
    DATE_FORMAT = '%b_%Y'
    TEXT = '.txt'
    WRITE_MODE = 'w'
    APPEND_MODE = 'a'
    SEPARATOR = '\t'
    NA = 'NA'
    EMPTY_STRING = ''


class FileDirectories():
    DATABASE_PATH = os.path.dirname(__file__)
    TEMP_DIR = 'temporary_outputs'
    UID_OUTPUT_DIR = 'temporary_outputs/entrez_id_list'
    ENTREZ_METADATA_DIR = 'temporary_outputs/entrez_species_metadata'
    ENTREZ_LOG_DIR = 'temporary_outputs/entrez_log'
    SPECIES_METADATA_DIR = 'temporary_outputs/species_reorganized_metadata'
    GC_METADATA_DIR = 'temporary_outputs/species_reorganized_metadata_gc'
    ASSEMBLY_DOWNLOAD_DIR = 'temporary_outputs/assembly_downloads'
    SPECIES_GC_LOG_DIR = 'temporary_outputs/species_gc_log'


class LogFiles():
    MAIN_LOG = 'LOG'
    SUB_LOG_ENTREZ = 'log_entrez_metadata_chunk'
    SUB_LOG_GC = '_log_gc'


class Assembly():
    ASSEMBLY_DATABASE = 'Assembly'
    CONTIG_TERM = 'contig[Assembly Level]'
    ORGANISM_AND_CONTIG_TERM = '[Organism] AND contig[Assembly Level]'
    RECORD_COUNT = 'Count'
    RECORD_IDLIST = 'IdList'
    ID_LIST_JOIN_CHAR = ','
    ESUMMARY_BATCH_LIMIT = 10000
    ESUMMARY_VALIDATE = False
    DOCUMENT_SUMMARY_SET = 'DocumentSummarySet'
    DOCUMENT_SUMMARY = 'DocumentSummary'
    ESUMMARY_UID_KEY = 'uid'
    SPECIES = 'SpeciesName'
    FTP_PATH_GENBANK = 'FtpPath_GenBank'
    ASSEMBLY_FILE_EXTENSION = '_genomic.fna.gz'
    ASSEMBLY_COUNTS_OUTPUT_PREFIX = 'species_assemblycounts_'
    ASSEMBLY_COUNT_LOWERBOUND = 10
    UID_PREFIX = 'Assembly_UID_chunk'
    METADATA_SUFFIX = '_metadata'
    COL_SPECIES = 'Species'
    COL_NUM_ASSEMBLIES = 'Num_assemblies'
    COL_NUM_UIDs = 'Num_UIDs'
    UID_SUFFIX = '_UIDs_numbers'
    EXCLUDED_SPECIES = 'excluded_species'
    EXCLUDED_ASSEMBLIES = 'excluded_assemblies'
    GC_SUFFIX = '_gc'


class Taxonomy():
    TAXONOMY_DATABASE = 'Taxonomy'
    TAXONOMY_RETMODE = 'xml'
    TAXONOMY_FILE_PREFIX = 'species_assemblycounts_'
    TAXONOMY_FILE_SUFFIX = '_taxonomy'
    LINEAGE = 'LineageEx'
    RANK = 'Rank'
    SCIENTIFIC_NAME = 'ScientificName'
    SUPERKINGDOM = 'superkingdom'
    SUPERKINGDOM_INDEX = 0
    PHYLUM = 'phylum'
    PHYLUM_INDEX = 1
    CLASS = 'class'
    CLASS_INDEX = 2
    ORDER = 'order'
    ORDER_INDEX = 3
    FAMILY = 'family'
    FAMILY_INDEX = 4
    GENUS = 'genus'
    GENUS_INDEX = 5
    PROKARYOTES = ['Bacteria', 'Archaea']
    TAXONOMY_WRITE_COLUMNS = ['Species', 'Num_assemblies', 'Kingdom', 'Phylum', 'Class', \
		'Order', 'Family', 'Genus']


class Metadata():
    METADATA_COLUMN_HEADERS = ['Organism Name', 
                            'Strain/Isolate', 
                            'Assembly Name', 
                            'Genbank Accession',
                            'Refseq Accession', 
                            'Genome Coverage',
                            'Submission Date', 
                            'Last Update Date',
                            'Refseq Exclusion Reason', 
                            'ContigN50', 
                            'Contig count', 
                            'ContigL50',
                            'Total length', 
                            'Assembly Method', 
                            'Sequencing Technology',
                            'GCcontent'
                            ]

    METADATA_INDEX_SPECIES = 0
    METADATA_INDEX_STRAIN = 1
    METADATA_INDEX_ASSEMBLY_ID = 2
    METADATA_INDEX_GENBANK = 3
    METADATA_INDEX_REFSEQ = 4
    METADATA_INDEX_COVERAGE = 5
    METADATA_INDEX_SUBMISSION = 6
    METADATA_INDEX_LAST_UPDATE = 7
    METADATA_INDEX_REFSEQ_EXCLUSION = 8
    METADATA_INDEX_N50 = 9
    METADATA_INDEX_NUM_CONTIGS = 10
    METADATA_INDEX_L50 = 11
    METADATA_INDEX_LENGTH = 12
    METADATA_INDEX_ASSEMBLER = 13
    METADATA_INDEX_SEQUENCING_TECHNOLOGY = 14
    METADATA_INDEX_GC_CONTENT = 15

