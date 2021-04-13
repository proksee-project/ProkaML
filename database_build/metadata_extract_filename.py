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
import sys
from Bio import Entrez
import urllib.request
import urllib.error
import re


#outfile.write('Organism Name\tOrganism Detail\tStrain/Isolate\tAssembly Name\tGenbank Accession\tRefseq Accession\tGenome Coverage\tSubmission Date\tLast Update Date\tRefseq Exclusion Reason\tContigN50\tContig count\tContigL50\tTotal length\tAssembly Method\tSequencing Technology\n')
def esummary(idlist=None, outfile=None):
	esum = Entrez.esummary(db="assembly", id=",".join(idlist))
	docsum = Entrez.read(esum, validate=False)
	for j in range(0, len(idlist)):
		d = docsum['DocumentSummarySet']['DocumentSummary'][j]
		org = d['SpeciesName']

		if 'Organism' in d:
			org_detail = d['Organism']
		else:
			org_detail = 'NA'

		if 'Biosource' in d:
			if d['Biosource']['InfraspeciesList']:
				subtype = d['Biosource']['InfraspeciesList'][0]['Sub_type']
				subval = d['Biosource']['InfraspeciesList'][0]['Sub_value']
				infra = subtype + ': ' + subval
			elif d['Biosource']['Isolate'] != '':
				isolate = d['Biosource']['Isolate']
				infra = 'isolate: ' + isolate
			else:
				infra = 'NA'
		else:
			infra = 'NA'

		if 'AssemblyName' in d:
			assembly_name = d['AssemblyName']
		else:
			assembly_name = 'NA'

		if 'Synonym' in d:
			if d['Synonym']['Genbank'] == '':
				genbank = 'NA'
			else:
				genbank = d['Synonym']['Genbank']
			if d['Synonym']['RefSeq'] == '':
				refseq = 'NA'
			else:
				refseq = d['Synonym']['RefSeq']
		else:
			genbank = 'NA'
			refseq = 'NA'

		if 'Coverage' in d:
			coverage = d['Coverage']
		else:
			coverage = 'NA'

		if 'SubmissionDate' in d:
			sub_date = d['SubmissionDate']
		else:
			sub_date = 'NA'

		if 'LastUpdateDate' in d:
			up_date = d['LastUpdateDate']
		else:
			up_date = 'NA'

		if 'ExclFromRefSeq' in d:
			refseq_excl = d['ExclFromRefSeq']
		else:
			refseq_excl = 'NA'

		if 'ContigN50' in d:
			n50 = d['ContigN50']
		else:
			n50 = 'NA'

		if 'Meta' in d:
			meta1 =  re.search(r'contig_count.+?(\d+)<', d['Meta'])
			if meta1 is not None:
				contig_count = meta1.group(1)
			else:
				contig_count = 'NA'
			meta2 =  re.search(r'contig_l50.+?(\d+)<', d['Meta'])
			if meta2 is not None:
				l50 = meta2.group(1)
			else:
				l50 = 'NA'
			meta3 =  re.search(r'total_length.+?(\d+)<', d['Meta'])
			if meta3 is not None:
				total_len = meta3.group(1)
			else:
				total_len = 'NA'
		else:
			contig_count = 'NA'
			l50 = 'NA'
			total_len = 'NA'

		if 'FtpPath_Assembly_rpt' in d and d['FtpPath_Assembly_rpt'] != '':
			url = d['FtpPath_Assembly_rpt']
			try:
				url_req = urllib.request.urlopen(url, timeout=100)

				# read() generates a string
				# decode() is applicable on a string
				# splitlines() removes \r and \n characters, generates a list
				lines = url_req.read().decode('utf-8').splitlines()
				lines_string = ';'.join(lines)
				meth = re.search(r'Assembly\smethod:\s(.+?);', lines_string)
				plat = re.search(r'Sequencing\stechnology:\s(.+?);', lines_string)

				if meth is not None:
					method = meth.group(1)
				else:
					method = 'NA'
				if plat is not None:
					platform = plat.group(1)
				else:
					platform = 'NA'

			except urllib.error.HTTPError:
				print('HTTP error. Ignoring')
			except urllib.error.URLError:
				print('URL error. NCBI update issue. Autofixes. Ignoring')

		else:
			method = 'NA'
			platform = 'NA'

		output_string = org + '\t' + org_detail + '\t' + infra + '\t' + assembly_name + '\t' + genbank + '\t' + refseq + '\t' + coverage + '\t' + sub_date + '\t' + up_date + '\t' + str(refseq_excl) + '\t' + n50 + '\t' + contig_count + '\t' + l50 + '\t' + total_len + '\t' + method + '\t' + platform
		outfile.write(output_string + '\n')
		print(str(j) + 'th record processed. GbUid: ' + d['GbUid'])


def idlist_parser(in_dir=None, out_dir=None, filename=None):
	with open(os.path.join(in_dir, filename)) as f:
		idlist = f.read().splitlines()
	
	output_filename = filename.split('idlist.txt')[0] + 'metadata.txt'
	output_file = open(os.path.join(out_dir, output_filename), 'w')
	
	return idlist, output_file


def main():
	if len(sys.argv) != 6:
		sys.exit('''
		Command usage: python assemblydb_entrez_query.py EMAIL NCBI_API_KEY INPUT_DIRECTORY OUTPUT_DIRECTORY FILE_NAME. 
		Need to pass 5 arguments corresponding to your email, ncbi api key, input directory of files containing
		UIDs, output dirctory and a file name containing list of UIDs. 
		''')

	else:
		email = sys.argv[1]
		api_key = sys.argv[2]
		input_dir = sys.argv[3]
		output_dir = sys.argv[4]
		filename = sys.argv[5]

		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

		Entrez.email = email
		Entrez.api_key = api_key
		idlist, output_file = idlist_parser(input_dir, output_dir, filename)
		esummary(idlist, output_file)


if __name__ == "__main__":
	main()