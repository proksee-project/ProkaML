from Bio import Entrez
from collections import defaultdict
import sys
import os


def count_assem_records(email=None, api_key=None):

	'''Declaring email and api key for making NCBI Entrez API calls'''
	Entrez.email = email
	Entrez.api_key = api_key

	try:
		'''Retrieving count of all contig assemblies'''
		handle = Entrez.esearch(db="assembly", term="contig[Assembly Level]")
		record = Entrez.read(handle)
		num = record['Count'] # num is approx 620k as of Sept 2020

	except Exception:
		raise Exception('Either internet failure OR invalid/incorrect ordering \
			of input parameters OR some unexpected failure..Exiting...')

		'''Returning total number of records'''
	return num


def retr_assem_records(num=None):

	'''Retrieving all records of NCBI contig assemblies'''
	handle1 = Entrez.esearch(db="assembly", term="contig[Assembly Level]", retmax = num)
	record1 = Entrez.read(handle1)

	'''esummary function can process maximum of 10,000 records in a single call
	So esummary is executed in batches of 10,000 records'''

	species_dicn = defaultdict(int)
	retstart = 0

	while retstart < int(num):
		idlist_batch = []
		for i in range(0, 10000):
			try:
				idlist_batch.append(record1['IdList'][i + retstart])
			except IndexError:
				last_batch_size = i
				break

		esummary(idlist_batch, species_dicn)

		try:
			retstart += last_batch_size
		except NameError:
			retstart += 10000

		str1 = str(retstart) + ' document summaries processed. '
		str2 = 'Species dictionary has ' + str(len(species_dicn)) + ' records'
		print (str1 + str2)

	return species_dicn


def esummary(idlist_batch=None, species_dicn=None):
	try:
		esum = Entrez.esummary(db="assembly", id=",".join(idlist_batch))
		docsum = Entrez.read(esum, validate=False)
		for j in range(0, len(idlist_batch)):
			species = docsum['DocumentSummarySet']['DocumentSummary'][j]['SpeciesName']
			species_dicn[species] += 1
	except Exception:
		raise Exception('Internet not working OR cannot parse esummary xml record')


def species_dicn_write(species_dicn=None, out_dir=None):
	output_file = open(os.path.join(out_dir, 'species_counts.txt'), 'w')

	for species in sorted(species_dicn, key=species_dicn.get, reverse=True):
		output_file.write('{}\t{}\n'.format(species, species_dicn[species]))


def main():
	email = sys.argv[1]
	api_key = sys.argv[2]
	output_dir = sys.argv[3]

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	num = count_assem_records(email, api_key)
	species_dicn = retr_assem_records(num)
	species_dicn_write(species_dicn, output_dir)

if __name__ == '__main__':
	main()