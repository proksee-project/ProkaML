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
import urllib.request
import urllib.error
import re
import constants as const


class EntrezMetadata():
    """
    A class for obtaining genomic assembly attributes from NCBI

    ATTRIBUTES
        idlist (list): list of NCBI assembly UIDs (int)
        document_summary (dict): nested biopython dictionary mapping assembly UIDs (int)
        to genomic attributes (str, float, int)
    """

    def __init__(self, idlist):
        """
        Initializes class for obtaining genomic metadata from NCBI contig assemblies

        PARAMETERS:
            idlist (list): list of NCBI assembly UIDs (int)
        """

        self.idlist = idlist

    def get_document_summaries(self, idlist):

        for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
            try:
                esum = Entrez.esummary(db=const.ASSEMBLY_DATABASE, id=const.ID_LIST_JOIN_CHAR.join(self.idlist))
                document_summaries = Entrez.read(esum, validate=const.ESUMMARY_VALIDATE)

            except Exception:
                document_summaries = {}

            else:
                if document_summaries:
                    break

        return document_summaries

    def print_genomic_metadata(self, outfile):
        """
        Prints/writes assembly metadata to output file

        PARAMETERS:
            outfile : the output file to which metadata is written

        POST
            Metadata for every assembly UID is written as column separated values (str, float, int)
            to output file
        """

        document_summaries = self.get_document_summaries(self.idlist)
        idlist_success_count = 0
        idlist_failures_count = 0
        irretrievable_uids = []
        successful_uids = []

        if len(document_summaries) > 0:
            idlist_success_count = len(document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY])
            if idlist_success_count == len(self.idlist):
                for j in range(0, len(self.idlist)):
                    document_dict = document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY][j]
                    metadata_string = const.SEPARATOR.join(self.get_metadata(document_dict))
                    outfile.write(metadata_string + const.NEW_LINE)
                    print('Count ' + str(int(j+1)) + ' assembly ' + document_dict['Synonym']['Genbank'] + ' : metadata obtained')

            else:
                idlist_failures_count = len(self.idlist) - len(document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY])
                for j in range(0, len(document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY])):
                    retrievable_uid = document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY][j].attributes['uid']
                    successful_uids.append(retrievable_uid)
                    document_dict = document_summaries[const.DOCUMENT_SUMMARY_SET][const.DOCUMENT_SUMMARY][j]
                    metadata_string = const.SEPARATOR.join(self.get_metadata(document_dict))
                    outfile.write(metadata_string + const.NEW_LINE)
                    print('Count ' + str(int(j+1)) + ' assembly ' + document_dict['Synonym']['Genbank'] + ' : metadata obtained')

                for k in range(0, len(self.idlist)):
                    if self.idlist[k] not in successful_uids:
                        irretrievable_uids.append(self.idlist[k])

        else:
            idlist_failures_count = len(self.idlist)
            irretrievable_uids.append(self.idlist)

        return idlist_success_count, idlist_failures_count, irretrievable_uids

    def get_metadata(self, document_dict):
        """
        Obtains metadata attributes for every assembly

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS
            metadata (list): list of assembly specific attributes (str, float, int)
        """

        species = self.get_species_name(document_dict)
        strain = self.get_species_strain(document_dict)
        assembly_id = self.get_assembly_id(document_dict)
        genbank_id = self.get_genbank_id(document_dict)
        refseq_id = self.get_refseq_id(document_dict)
        genome_coverage = self.get_genome_coverage(document_dict)
        submission_date = self.get_submission_date(document_dict)
        last_update_date = self.get_last_update_date(document_dict)
        refseq_exclusion_reason = str(self.get_refseq_exclusion_reason(document_dict))
        n50 = self.get_n50(document_dict)
        num_contigs = self.get_num_contigs(document_dict)
        l50 = self.get_l50(document_dict)
        length = self.get_length(document_dict)
        assembly_report = self.get_assembly_report(document_dict)

        assembler = self.get_assembler(assembly_report)
        sequencing_platform = self.get_sequencing_platform(assembly_report)

        metadata = [species, strain, assembly_id, genbank_id, refseq_id, genome_coverage, submission_date,
                    last_update_date, refseq_exclusion_reason, n50, num_contigs, l50, length, assembler,
                    sequencing_platform]

        return metadata

    def get_species_name(self, document_dict):
        """
        Obtains species name for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            species_name (str): the species name for an assembly from NCBI
        """

        if 'SpeciesName' in document_dict:
            species_name = document_dict['SpeciesName']

        else:
            species_name = 'NA'

        return species_name

    def get_species_strain(self, document_dict):
        """
        Obtains species strain/isolate for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            strain (str): the species strain/isolate (if applicable) for an assembly from NCBI
        """

        if 'Biosource' in document_dict:
            if document_dict['Biosource']['InfraspeciesList']:
                subtype = document_dict['Biosource']['InfraspeciesList'][0]['Sub_type']
                subvalue = document_dict['Biosource']['InfraspeciesList'][0]['Sub_value']
                strain = subtype + ': ' + subvalue

            elif document_dict['Biosource']['Isolate'] != '':
                isolate = document_dict['Biosource']['Isolate']
                strain = 'isolate: ' + isolate

            else:
                strain = 'NA'

        else:
            strain = 'NA'

        return strain

    def get_assembly_id(self, document_dict):
        """
        Obtains assembly ID for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            assembly_id (str): the assembly ID for an assembly from NCBI
        """

        if 'AssemblyName' in document_dict:
            assembly_id = document_dict['AssemblyName']

        else:
            assembly_id = 'NA'

        return assembly_id

    def get_genbank_id(self, document_dict):
        """
        Obtains Genbank ID for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            genbank_id (str): the Genbank ID for an assembly from NCBI
        """

        if 'Synonym' in document_dict:
            if document_dict['Synonym']['Genbank'] == '':
                genbank_id = 'NA'

            else:
                genbank_id = document_dict['Synonym']['Genbank']

        else:
            genbank_id = 'NA'

        return genbank_id

    def get_refseq_id(self, document_dict):
        """
        Obtains RefSeq Accession ID for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            refseq_id (str): the RefSeq Accession ID (if applicable) for an assembly from NCBI
        """

        if 'Synonym' in document_dict:
            if document_dict['Synonym']['RefSeq'] == '':
                refseq_id = 'NA'

            else:
                refseq_id = document_dict['Synonym']['RefSeq']

        else:
            refseq_id = 'NA'

        return refseq_id

    def get_genome_coverage(self, document_dict):
        """
        Obtains genomic coverage for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            genome_coverage (float): the genomic coverage for an assembly from NCBI
        """

        if 'Coverage' in document_dict:
            genome_coverage = document_dict['Coverage']

        else:
            genome_coverage = 'NA'

        return genome_coverage

    def get_submission_date(self, document_dict):
        """
        Obtains submission date for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            submission_date (str): the submission date for an assembly from NCBI
        """

        if 'SubmissionDate' in document_dict:
            submission_date = document_dict['SubmissionDate']

        else:
            submission_date = 'NA'

        return submission_date

    def get_last_update_date(self, document_dict):
        """
        Obtains last updated date for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            last_update_date (str): the last update date for an assembly from NCBI
        """

        if 'LastUpdateDate' in document_dict:
            last_update_date = document_dict['LastUpdateDate']

        else:
            last_update_date = 'NA'

        return last_update_date

    def get_refseq_exclusion_reason(self, document_dict):
        """
        Obtains RefSeq exclusion reason/s (if applicable) for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            refseq_excl (str): RefSeq exclusion reason/s (if applicable) for an assembly from NCBI
        """

        if 'ExclFromRefSeq' in document_dict:
            refseq_excl = document_dict['ExclFromRefSeq']

        else:
            refseq_excl = 'NA'

        return refseq_excl

    def get_n50(self, document_dict):
        """
        Obtains N50 for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            N50 (int): the N50 for an assembly from NCBI
        """

        if 'ContigN50' in document_dict:
            n50 = document_dict['ContigN50']

        else:
            n50 = float('NaN')

        return n50

    def get_num_contigs(self, document_dict):
        """
        Obtains number of contigs for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            num_contigs (int): the number of contigs for an assembly from NCBI
        """

        if 'Meta' in document_dict:
            num_contigs_regex = re.search(r'contig_count.+?(\d+)<', document_dict['Meta'])
            if num_contigs_regex is not None:
                num_contigs = num_contigs_regex.group(1)

            else:
                num_contigs = float('NaN')

        else:
            num_contigs = float('NaN')

        return num_contigs

    def get_l50(self, document_dict):
        """
        Obtains L50 for an assembly from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            l50 (int): the L50 for an assembly from NCBI
        """

        if 'Meta' in document_dict:
            l50_regex = re.search(r'contig_l50.+?(\d+)<', document_dict['Meta'])
            if l50_regex is not None:
                l50 = l50_regex.group(1)

            else:
                l50 = float('NaN')

        else:
            l50 = float('NaN')

        return l50

    def get_length(self, document_dict):
        """
        Obtains assembly length from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            length (int): the assembly length from NCBI
        """

        if 'Meta' in document_dict:
            length_regex = re.search(r'total_length.+?(\d+)<', document_dict['Meta'])
            if length_regex is not None:
                length = length_regex.group(1)

            else:
                length = float('NaN')

        else:
            length = float('NaN')

        return length

    def get_assembly_report(self, document_dict):
        """
        Obtains assembly report length from NCBI

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS:
            assembly_report_lines_string (str): the assembly report from NCBI
        """

        ASSEMBLY_REPORT_KEY = 'FtpPath_Assembly_rpt'
        TIMEOUT = 100
        JOIN_CHARACTER = ';'

        if ASSEMBLY_REPORT_KEY in document_dict and document_dict[ASSEMBLY_REPORT_KEY] != '':
            assembly_report_url = document_dict[ASSEMBLY_REPORT_KEY]

            try:
                assembly_url_request = urllib.request.urlopen(assembly_report_url, timeout=TIMEOUT)

                # read() generates a string
                # decode() is applicable on a string
                # splitlines() removes \r and \n characters, generates a list
                assembly_report_lines_list = assembly_url_request.read().decode('utf-8').splitlines()
                assembly_report_lines_string = JOIN_CHARACTER.join(assembly_report_lines_list)

            except (urllib.error.HTTPError, urllib.error.URLError):
                """
                url errors can arise due to NCBI server updates. The timeout argument in the try: block
                should accomodate most autofixing of url error issues.
                """
                assembly_report_lines_string = 'URLError'

        else:
            assembly_report_lines_string = 'NA'

        return assembly_report_lines_string

    def get_assembler(self, assembly_report_lines_string):
        """
        Obtains assembly method from NCBI

        PARAMETERS:
            assembly_report_lines_string (str): the assembly report from NCBI

        RETURNS:
            assembler (str): the assembly method from NCBI
        """

        if assembly_report_lines_string == 'NA' or assembly_report_lines_string == 'URLError':
            assembler = 'NA'

        else:
            assembler_regex = re.search(r'Assembly\smethod:\s(.+?);', assembly_report_lines_string)
            if assembler_regex is not None:
                assembler = assembler_regex.group(1)

            else:
                assembler = 'NA'

        return assembler

    def get_sequencing_platform(self, assembly_report_lines_string):
        """
        Obtains sequencing platform from NCBI

        PARAMETERS:
            assembly_report_lines_string (str): the assembly report from NCBI

        RETURNS:
            sequencing_platform (str): the sequencing platform for an assembly from NCBI
        """

        if assembly_report_lines_string == 'NA' or assembly_report_lines_string == 'URLError':
            sequencing_platform = 'NA'
 
        else:
            sequencing_platform_regex = re.search(r'Sequencing\stechnology:\s(.+?);', assembly_report_lines_string)
            if sequencing_platform_regex is not None:
                sequencing_platform = sequencing_platform_regex.group(1)

            else:
                sequencing_platform = 'NA'

        return sequencing_platform