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

from Bio import Entrez
import urllib.request
import urllib.error
import re


class GetMetadata():
    """
    A class for obtaining genomic assembly attributes from NCBI

    ATTRIBUTES
        idlist (list): list of NCBI assembly UIDs
    """

    def __init__(self, idlist):
        """
        Initializes class for obtaining genomic metadata from NCBI contig assemblies

        PARAMETERS:
            idlist (list): list of NCBI assembly UIDs
        """

        self.idlist = idlist

        # Obtain document summaries of assemblies using Entrez esummary function
        esum = Entrez.esummary(db="assembly", id=",".join(self.idlist))
        self.document_summary = Entrez.read(esum, validate=False)

    def print_genomic_metadata(self, outfile):
        """
        Prints/writes assembly metadata to output file

        PARAMETERS:
            outfile : the output file to which metadata is written

        POST
            Metadata for every assembly UID is written as column separated entities to output file
        """

        for j in range(0, len(self.idlist)):
            document_dict = self.document_summary['DocumentSummarySet']['DocumentSummary'][j]

            separator = '\t'
            metadata_string = separator.join(self.get_metadata(document_dict))
            outfile.write(metadata_string + '\n')
            print(str(j) + 'th record processed. GbUid: ' + document_dict['GbUid'])

    def get_metadata(self, document_dict):
        """
        Obtains metadata attributes for every assembly

        PARAMETERS:
            document_dict (dict): the document dictionary of an assembly obtained by NCBI Entrez
            esummary function
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

        if 'SpeciesName' in document_dict:
            species_name = document_dict['SpeciesName']

        else:
            species_name = 'NA'

        return species_name

    def get_species_strain(self, document_dict):

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

        if 'AssemblyName' in document_dict:
            assembly_id = document_dict['AssemblyName']

        else:
            assembly_id = 'NA'

        return assembly_id

    def get_genbank_id(self, document_dict):

        if 'Synonym' in document_dict:
            if document_dict['Synonym']['Genbank'] == '':
                genbank_id = 'NA'

            else:
                genbank_id = document_dict['Synonym']['Genbank']

        else:
            genbank_id = 'NA'

        return genbank_id

    def get_refseq_id(self, document_dict):

        if 'Synonym' in document_dict:
            if document_dict['Synonym']['RefSeq'] == '':
                refseq_id = 'NA'

            else:
                refseq_id = document_dict['Synonym']['RefSeq']

        else:
            refseq_id = 'NA'

        return refseq_id

    def get_genome_coverage(self, document_dict):

        if 'Coverage' in document_dict:
            genome_coverage = document_dict['Coverage']

        else:
            genome_coverage = 'NA'

        return genome_coverage

    def get_submission_date(self, document_dict):

        if 'SubmissionDate' in document_dict:
            submission_date = document_dict['SubmissionDate']

        else:
            submission_date = 'NA'

        return submission_date

    def get_last_update_date(self, document_dict):

        if 'LastUpdateDate' in document_dict:
            last_update_date = document_dict['LastUpdateDate']

        else:
            last_update_date = 'NA'

        return last_update_date

    def get_refseq_exclusion_reason(self, document_dict):

        if 'ExclFromRefSeq' in document_dict:
            refseq_excl = document_dict['ExclFromRefSeq']

        else:
            refseq_excl = 'NA'

        return refseq_excl

    def get_n50(self, document_dict):

        if 'ContigN50' in document_dict:
            n50 = document_dict['ContigN50']

        else:
            n50 = 'NA'

        return n50

    def get_num_contigs(self, document_dict):

        if 'Meta' in document_dict:
            num_contigs_regex = re.search(r'contig_count.+?(\d+)<', document_dict['Meta'])
            if num_contigs_regex is not None:
                num_contigs = num_contigs_regex.group(1)

            else:
                num_contigs = 'NA'

        else:
            num_contigs = 'NA'

        return num_contigs

    def get_l50(self, document_dict):

        if 'Meta' in document_dict:
            l50_regex = re.search(r'contig_l50.+?(\d+)<', document_dict['Meta'])
            if l50_regex is not None:
                l50 = l50_regex.group(1)

            else:
                l50 = 'NA'

        else:
            l50 = 'NA'

        return l50

    def get_length(self, document_dict):

        if 'Meta' in document_dict:
            length_regex = re.search(r'total_length.+?(\d+)<', document_dict['Meta'])
            if length_regex is not None:
                length = length_regex.group(1)

            else:
                length = 'NA'

        else:
            length = 'NA'

        return length

    def get_assembly_report(self, document_dict):

        if 'FtpPath_Assembly_rpt' in document_dict and document_dict['FtpPath_Assembly_rpt'] != '':
            assembly_report_url = document_dict['FtpPath_Assembly_rpt']

            try:
                assembly_url_request = urllib.request.urlopen(assembly_report_url, timeout=100)

                # read() generates a string
                # decode() is applicable on a string
                # splitlines() removes \r and \n characters, generates a list
                assembly_report_lines_list = assembly_url_request.read().decode('utf-8').splitlines()
                assembly_report_lines_string = ';'.join(assembly_report_lines_list)

            except (urllib.error.HTTPError, urllib.error.URLError):
                """
                url errors can arise due to NCBI server updates. The timeout argument in the try: block
                is sufficiently large to accommodate autofixing of url error issues.
                """

                pass

        else:
            assembly_report_lines_string = 'NA'

        return assembly_report_lines_string

    def get_assembler(self, assembly_report_lines_string):

        if assembly_report_lines_string == 'NA':
            assembler = 'NA'

        else:
            assembler_regex = re.search(r'Assembly\smethod:\s(.+?);', assembly_report_lines_string)
            if assembler_regex is not None:
                assembler = assembler_regex.group(1)

            else:
                assembler = 'NA'

        return assembler

    def get_sequencing_platform(self, assembly_report_lines_string):

        if assembly_report_lines_string == 'NA':
            sequencing_platform = 'NA'

        else:
            sequencing_platform_regex = re.search(r'Sequencing\stechnology:\s(.+?);', assembly_report_lines_string)
            if sequencing_platform_regex is not None:
                sequencing_platform = sequencing_platform_regex.group(1)

            else:
                sequencing_platform = 'NA'

        return sequencing_platform
