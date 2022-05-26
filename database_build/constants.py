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

Entrez.max_tries = 1
Entrez.sleep_between_tries = 1
API_QUERY_ATTEMPT_START = 1
API_QUERY_ATTEMPT_END = 4

DATE_FORMAT = '%b_%Y'
LOG_FILE = 'LOG'
FILE_EXTENSION = '.txt'
WRITE_MODE = 'w'
APPEND_MODE = 'a'
SEPARATOR = '\t'
NEW_LINE = '\n'
NA = 'NA'
UID_PREFIX = 'Assembly_UID_chunk'
METADATA_SUFFIX = '_metadata'
EMPTY_STRING = ''

COL_KINGDOM = 'Kingdom'
COL_SPECIES = 'Species'
COL_NUM_ASSEMBLIES = 'Num_assemblies'
COL_NUM_UIDs = 'Num_UIDs'

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
ASSEMBLY_COUNTS_OUTPUT_PREFIX = 'species_assemblycounts_'
ASSEMBLY_COUNTS_OUTPUT_HEADER = 'Species\tNum_assemblies\n'
ASSEMBLY_COUNT_LOWERBOUND = 10

TAXONOMY_DATABASE = 'Taxonomy'
TAXONOMY_RETMODE = 'xml'
TAXONOMY_FILE_PREFIX = ASSEMBLY_COUNTS_OUTPUT_PREFIX
TAXONOMY_FILE_SUFFIX = '_taxonomy'
LINEAGE = 'LineageEx'
RANK = 'Rank'
SCIENTIFIC_NAME = 'ScientificName'
SUPERKINGDOM = 'superkingdom'
PHYLUM = 'phylum'
CLASS = 'class'
ORDER = 'order'
FAMILY = 'family'
GENUS = 'genus'

UID_SUFFIX = '_UIDs_numbers'
PROKARYOTES = ['Bacteria', 'Archaea']
UID_OUTPUT_DIR = 'entrez_id_list'
ENTREZ_METADATA_DIR = 'entrez_species_metadata'
REORGANIZED_METADATA_DIR = 'species_reorganized_metadata'
ADDITIONAL_METADATA_DIR = 'additional_species_metadata'

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
                           'Sequencing Technology'
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
