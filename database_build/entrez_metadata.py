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


class EntrezMetadataDownloader():
    """
    A class for obtaining genomic assembly attributes from NCBI

    ATTRIBUTES
        idlist (list): list of NCBI assembly UIDs (int)
        output_file (str): assembly chunk specific output file
    """

    def __init__(self, idlist, output_file):
        """
        Initializes class for obtaining genomic metadata from NCBI contig assemblies

        PARAMETERS:
            idlist (list): list of NCBI assembly UIDs (int)
            output_file (str): assembly chunk specific output file
        """

        self.idlist = idlist
        self.output_file = output_file

    def __get_document_summaries(self, idlist):
        """
        Obtains NCBI document summaries for a list of assembly UIDs

        PARAMETERS:
            idlist (list): list of NCBI assembly UIDs (int)
        
        RETURNS:
            document_summaries (dict): nested dictionary of assembly attributes from NCBI
        """

        for attempts in range(const.API_QUERY_ATTEMPT_START, const.API_QUERY_ATTEMPT_END):
            try:
                esum = Entrez.esummary(db=const.Assembly.ASSEMBLY_DATABASE, id=const.Assembly.ID_LIST_JOIN_CHAR.join(self.idlist))
                document_summaries = Entrez.read(esum, validate=const.Assembly.ESUMMARY_VALIDATE)

            except Exception:
                document_summaries = {}

            else:
                if document_summaries:
                    break

        return document_summaries

    def __get_metadata(self, document_dict):
        """
        Obtains metadata attributes for every assembly

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly mapping UID (int)
            to assembly attributes (str, float, int)

        RETURNS
            metadata (list): list of assembly specific attributes (str, float, int)
        """

        species = self.__get_species_name(document_dict)
        strain = self.__get_species_strain(document_dict)
        assembly_id = self.__get_assembly_id(document_dict)
        genbank_id = self.__get_genbank_id(document_dict)
        refseq_id = self.__get_refseq_id(document_dict)
        genome_coverage = self.__get_genome_coverage(document_dict)
        submission_date = self.__get_submission_date(document_dict)
        last_update_date = self.__get_last_update_date(document_dict)
        refseq_exclusion_reason = str(self.__get_refseq_exclusion_reason(document_dict))
        n50 = self.__get_n50(document_dict)
        num_contigs = self.__get_num_contigs(document_dict)
        l50 = self.__get_l50(document_dict)
        length = self.__get_length(document_dict)
        assembly_report = self.__get_assembly_report(document_dict)

        assembler = self.__get_assembler(assembly_report)
        sequencing_platform = self.__get_sequencing_platform(assembly_report)

        metadata = [species, strain, assembly_id, genbank_id, refseq_id, genome_coverage, submission_date,
                    last_update_date, refseq_exclusion_reason, n50, num_contigs, l50, length, assembler,
                    sequencing_platform]

        return metadata

    def __get_species_name(self, document_dict):
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
            species_name = const.FileFormat.NA

        return species_name

    def __get_species_strain(self, document_dict):
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
                strain = const.FileFormat.NA

        else:
            strain = const.FileFormat.NA

        return strain

    def __get_assembly_id(self, document_dict):
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
            assembly_id = const.FileFormat.NA

        return assembly_id

    def __get_genbank_id(self, document_dict):
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
                genbank_id = const.FileFormat.NA

            else:
                genbank_id = document_dict['Synonym']['Genbank']

        else:
            genbank_id = const.FileFormat.NA

        return genbank_id

    def __get_refseq_id(self, document_dict):
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
                refseq_id = const.FileFormat.NA

            else:
                refseq_id = document_dict['Synonym']['RefSeq']

        else:
            refseq_id = const.FileFormat.NA

        return refseq_id

    def __get_genome_coverage(self, document_dict):
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
            genome_coverage = const.FileFormat.NA

        return genome_coverage

    def __get_submission_date(self, document_dict):
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
            submission_date = const.FileFormat.NA

        return submission_date

    def __get_last_update_date(self, document_dict):
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
            last_update_date = const.FileFormat.NA

        return last_update_date

    def __get_refseq_exclusion_reason(self, document_dict):
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
            refseq_excl = const.FileFormat.NA

        return refseq_excl

    def __get_n50(self, document_dict):
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

    def __get_num_contigs(self, document_dict):
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

    def __get_l50(self, document_dict):
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

    def __get_length(self, document_dict):
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

    def __get_assembly_report(self, document_dict):
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
                assembly_report_lines_list = assembly_url_request.read().decode('utf-8').splitlines()
                assembly_report_lines_string = JOIN_CHARACTER.join(assembly_report_lines_list)

            except (urllib.error.HTTPError, urllib.error.URLError):
                """
                url errors can arise due to NCBI server updates. The timeout argument in the try: block
                should accomodate most autofixing of url error issues.
                """
                assembly_report_lines_string = 'URLError'

        else:
            assembly_report_lines_string = const.FileFormat.NA

        return assembly_report_lines_string

    def __get_assembler(self, assembly_report_lines_string):
        """
        Obtains assembly method from NCBI

        PARAMETERS:
            assembly_report_lines_string (str): the assembly report from NCBI

        RETURNS:
            assembler (str): the assembly method from NCBI
        """

        if assembly_report_lines_string == const.FileFormat.NA or assembly_report_lines_string == 'URLError':
            assembler = const.FileFormat.NA

        else:
            assembler_regex = re.search(r'Assembly\smethod:\s(.+?);', assembly_report_lines_string)
            if assembler_regex is not None:
                assembler = assembler_regex.group(1)

            else:
                assembler = const.FileFormat.NA

        return assembler

    def __get_sequencing_platform(self, assembly_report_lines_string):
        """
        Obtains sequencing platform from NCBI

        PARAMETERS:
            assembly_report_lines_string (str): the assembly report from NCBI

        RETURNS:
            sequencing_platform (str): the sequencing platform for an assembly from NCBI
        """

        if assembly_report_lines_string == const.FileFormat.NA or assembly_report_lines_string == 'URLError':
            sequencing_platform = const.FileFormat.NA
 
        else:
            sequencing_platform_regex = re.search(r'Sequencing\stechnology:\s(.+?);', assembly_report_lines_string)
            if sequencing_platform_regex is not None:
                sequencing_platform = sequencing_platform_regex.group(1)

            else:
                sequencing_platform = const.FileFormat.NA

        return sequencing_platform

    def __print_genomic_metadata(self, document_summaries, document_summary_index):
        """
        Prints/writes assembly metadata to output file

        PARAMETERS:
            document_summaries (dict): nested dictionary of assembly attributes from NCBI
            document_summary_index (int): index of a document summary

        RETURNS:
            retrievable_uid (int): Assembly UID for which metadata is obtained
        """

        document_dict = document_summaries[const.Assembly.DOCUMENT_SUMMARY_SET][const.Assembly.DOCUMENT_SUMMARY][document_summary_index]
        metadata_string = const.FileFormat.SEPARATOR.join(self.__get_metadata(document_dict))
        self.output_file.write(metadata_string + '\n')
        print('Count ' + str(int(document_summary_index+1)) + ' assembly ' + document_dict['Synonym']['Genbank'] + \
            ' : metadata obtained')
        retrievable_uid = document_dict.attributes[const.Assembly.ESUMMARY_UID_KEY]

        return retrievable_uid

    def execute(self):
        """
        Prints/writes assembly metadata to output file

        RETURNS:
            idlist_success_count (int): number of UIDs for which metadata is successfully retrieved
            irretrievable_uids (list): list of UIDs for which metadata could not be retrieved
        """

        document_summaries = self.__get_document_summaries(self.idlist)
        idlist_success_count = 0
        irretrievable_uids = []
        successful_uids = []

        if len(document_summaries) > 0:
            idlist_success_count = len(document_summaries[const.Assembly.DOCUMENT_SUMMARY_SET][const.Assembly.DOCUMENT_SUMMARY])
            if idlist_success_count == len(self.idlist):
                for j in range(0, len(self.idlist)):
                    self.__print_genomic_metadata(document_summaries, j)

            else:
                for j in range(0, idlist_success_count):
                    retrievable_uid = self.__print_genomic_metadata(document_summaries, j)
                    successful_uids.append(retrievable_uid)

                for k in range(0, len(self.idlist)):
                    if self.idlist[k] not in successful_uids:
                        irretrievable_uids.append(self.idlist[k])

        else:
            irretrievable_uids.append(self.idlist)

        return idlist_success_count, irretrievable_uids
