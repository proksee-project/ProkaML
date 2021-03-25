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
from collections import defaultdict
import sys
import os
from datetime import date


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


def retrieve_assem_records(num=None):
	'''Retrieving all records of NCBI contig assemblies'''
	handle1 = Entrez.esearch(db="assembly", term="contig[Assembly Level]", retmax = num)
	record1 = Entrez.read(handle1)

	'''Initiating a species/organism counting dictionary'''
	species_dicn = defaultdict(int)

	'''esummary function can process maximum of 10,000 records in a single API call.
	So esummary is executed in batches of 10,000 records. 
	Batch incrementing is done by retstart variable'''
	retstart = 0

	while retstart < int(num):
		idlist_batch = []
		for i in range(0, 10000):
			try:
				'''10,000 UIDs are processed per batch'''
				idlist_batch.append(record1['IdList'][i + retstart])

				'''Catching exception for less than 10,000 UIDs for last batch'''
			except IndexError:
				last_batch_size = i
				break

				'''Calling esummary function'''
		esummary(idlist_batch, species_dicn)

		'''Incrementing batch count initiation by 10,000. Outer while loop exits for last batch'''
		try:
			retstart += last_batch_size
		except NameError:
			retstart += 10000

			'''Printing intermediate progresses''' 
		str1 = str(retstart) + ' document summaries processed. '
		str2 = 'Species dictionary has ' + str(len(species_dicn)) + ' records'
		print (str1 + str2)

		'''Returning dictionary of organism names and counts in assembly database'''
	return species_dicn


def esummary(idlist_batch=None, species_dicn=None):
	'''Calling Entrez esummary to generate xml document summary'''
	esum = Entrez.esummary(db="assembly", id=",".join(idlist_batch))

	'''Converting xml to python nested dictionary using Entrez read function'''
	docsum = Entrez.read(esum, validate=False)
	for j in range(0, len(idlist_batch)):

		'''Retrieving species name from nested dictionary'''
		species = docsum['DocumentSummarySet']['DocumentSummary'][j]['SpeciesName']
		species_dicn[species] += 1


def species_dicn_write(species_dicn=None):

	'''Writing species dictionary (sorted by descending occurences) to output file stamped with month and year'''
	month_year_stamp = date.today().strftime("%b_%Y")
	output_file = open('species_counts_' + month_year_stamp + '.txt', 'w')
	if species_dicn:
		for species in sorted(species_dicn, key=species_dicn.get, reverse=True):
			output_file.write('{}\t{}\n'.format(species, species_dicn[species]))
		success = 'Species dictionary written'
	else:
		success = 'No species dictionary from previous steps'
	
	return success


def main():
	if len(sys.argv) != 3:
		sys.exit('''
		Command usage: python assemblydb_entrez_query.py EMAIL NCBI_API_KEY. 
		Need to pass 2 arguments corresponding to your email and ncbi api key
		''')

	else:
		email = sys.argv[1]
		api_key = sys.argv[2]

		try:
			'''Obtaining total counts of assembly records'''
			num = count_assem_records(email, api_key)
		except Exception as e:
			sys.exit(e)

		''''Processing assembly records to create species dictionary'''
		species_dicn = retrieve_assem_records(num)

		'''Writing species dictionary to output'''
		print(species_dicn_write(species_dicn))


if __name__ == '__main__':
	main()
