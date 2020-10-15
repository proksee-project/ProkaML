from Bio import Entrez
import urllib.request
import urllib.error
import re

#outfile.write('Organism Name\tOrganism Detail\tAssembly Name\tGenbank Accession\tRefseq Accession\tGenome Coverage\tSubmission Date\tLast Update Date\tRefseq Exclusion Reason\tContigN50\tAssembly Method\tSequencing Technology\n')
class Metadata():

	def __init__(self, idlist, outfile):
		self.idlist = idlist
		self.outfile = outfile

	def docsum_dicn(self, d):
		org = d['SpeciesName']
		if 'Organism' in d:
			org_detail = d['Organism']
		else:
			org_detail = 'NA'
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

		output_string = org + '\t' + org_detail + '\t' + assembly_name + '\t' + genbank + '\t' + refseq + '\t' + coverage + '\t' + sub_date + '\t' + up_date + '\t' + str(refseq_excl) + '\t' + n50 + '\t'
		
		return output_string

	def docsum_meta(self, d):
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

		output_string = contig_count + '\t' + l50 + '\t' + total_len + '\t'

		return output_string

	def docsum_url(self, d):
		if 'FtpPath_Assembly_rpt' in d and d['FtpPath_Assembly_rpt'] != '':
			url = d['FtpPath_Assembly_rpt']
			try:
				url_req = urllib.request.urlopen(url, timeout=10)

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

		output_string = method + '\t' + platform
			
		return output_string

	def meta_integrate(self):
		esum = Entrez.esummary(db="assembly", id=",".join(self.idlist))
		docsum = Entrez.read(esum, validate=False)
		for j in range(0, len(self.idlist)):
			mono_doc = docsum['DocumentSummarySet']['DocumentSummary'][j]
			meta = self.docsum_dicn(mono_doc) + self.docsum_meta(mono_doc) + self.docsum_url(mono_doc)
			self.outfile.write(meta + '\n')
			print(str(j) + 'th record processed. GbUid: ' + mono_doc['GbUid'])
