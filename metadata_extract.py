import os
from Bio import Entrez
import urllib
import re

Entrez.email = "arnab22.iitkgp@gmail.com"
Entrez.api_key = "9b946150b1df235cf0e7c288ca4588a97c08"

query_sp = []
infile = open('species_counts.txt','r')
for line in infile:
	sp = line.rstrip().split('\t')
	if int(sp[1]) > 1000:
		query_sp.append(sp[0])

print(query_sp)
outfile = open('species_metadata.txt', 'w')
outfile.write('Organism Name\tAssembly Name\tGenbank Accession\tRefseq Accession\tGenome Coverage\tPartial Genome Represenation\tLast Update Date\tRefseq Exclusion Reason\tContigN50\tAssembly Method\tSequencing Technology\n')

for org in query_sp:
	handle = Entrez.esearch(db="assembly", term=org + "[Organism] AND contig[Assembly Level]", retmax=100)
	record = Entrez.read(handle)
	esum = Entrez.esummary(db="assembly", id=",".join(record['IdList']))
	docsum = Entrez.read(esum, validate=False)

	for i in range(0, len(record['IdList'])):
		d = docsum['DocumentSummarySet']['DocumentSummary'][i]

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

		if 'PartialGenomeRepresentation' in d:
			partial_genome_repr = d['PartialGenomeRepresentation']		
		else:
			partial_genome_repr = 'NA'

		if 'LastUpdateDate' in d:
			date = d['LastUpdateDate']
		else:
			date = 'NA'

		if 'ExclFromRefSeq' in d:
			refseq_excl = d['ExclFromRefSeq']
		else:
			refseq_excl = 'NA'

		if 'ContigN50' in d:
			n50 = d['ContigN50']
		else:
			n50 = 'NA'

		if 'FtpPath_Assembly_rpt' in d:
			url = d['FtpPath_Assembly_rpt']
			url_req = urllib.request.urlopen(url)

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

		else:
			method = 'NA'
			platform = 'NA'
		
		output_string = org + '\t' + assembly_name + '\t' + genbank + '\t' + refseq + '\t' + coverage + '\t' + partial_genome_repr + '\t' + date + '\t' + str(refseq_excl) + '\t' + n50 + '\t' + method + '\t' + platform + '\n'
		outfile.write(output_string)


	print(org + ' metadata written')