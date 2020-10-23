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
import re

def species_select(in_dir=None):
	sp_list_major = []
	sp_list_large = []
	sp_list_interm = []
	sp_list_minor = []
	
	infile = open(os.path.join(in_dir, 'species_counts_oct2020.txt'),'r')
	for line in infile:
		sp = line.rstrip().split('\t')
		sp_corrected = re.sub(r'[\[\]]','', sp[0])
		sp_tup = sp_corrected, sp[1]

		if int(sp_tup[1]) >= 1000:
			sp_list_major.append(sp_tup)
		elif int(sp_tup[1]) < 1000 and int(sp_tup[1]) >= 100:
			sp_list_large.append(sp_tup)
		elif int(sp_tup[1]) < 100 and int(sp_tup[1]) >= 10:
			sp_list_interm.append(sp_tup)
		elif int(sp_tup[1]) < 10:
			sp_list_minor.append(sp_tup)

	species_entire_tup = sp_list_major, sp_list_large, sp_list_interm, sp_list_minor

	return species_entire_tup

def id_list(species_tuple=None):
	try:
		handle = Entrez.esearch(db="assembly",
								term=species_tuple[0] + "[Organism] AND contig[Assembly Level]",
								retmax=species_tuple[1])
		record = Entrez.read(handle)

	except Exception:
		raise Exception('Either internet failure OR some unexpected failure..Exiting...')

	return record['IdList']

def chunk_counter(idlist=None):
	if len(idlist) % 10000 == 0:
		chunks = len(idlist)/10000
	else:
		chunks = int(len(idlist)/10000) + 1

	return chunks

def species_outfile(species_tuple=None, chunks=None):
	species_outfile_list = []
	species_name_list = re.split('\s|\/', species_tuple[0])
	for i in range(1, chunks + 1):
		outfile = "_".join(species_name_list) + '_chunk' + str(i) + '_idlist.txt'
		species_outfile_list.append(outfile)

	return species_outfile_list

def species_iterate(species_list=None, out_dir=None):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	for count, sp_tup in enumerate(species_list):
		idlist = id_list(sp_tup)
		chunks = chunk_counter(idlist)
		if chunks > 0:
			species_outfile_list = species_outfile(sp_tup, chunks)
			print(str(count) + ' species processing: ' + str(sp_tup[0]))
			for k in range(0, chunks):
				idlist_lowerlim = k*10000
				idlist_upperlim = (k+1)*10000
				idlist_batch = []
				for m in range(idlist_lowerlim, min(idlist_upperlim, len(idlist))):
						idlist_batch.append(idlist[m])
				output_handle = os.path.join(out_dir, species_outfile_list[k])
				output_file = open(output_handle, 'w')
				idlist_string = '\n'.join(idlist_batch)
				output_file.write(idlist_string +'\n')
		elif chunks == 0:
			print(str(count) + ' species processing: ' + str(sp_tup[0]) + ' returns empty ID list with esearch. skipping..')

def main():
	if len(sys.argv) != 4:
		sys.exit('''
		Command usage: python idlist_retriever.py EMAIL NCBI_API_KEY INPUT_DIRECTORY. 
		Need to pass 3 arguments corresponding to your email, your ncbi api key and 
		the input directory containing species_counts_timeofyear.txt
		''')
	
	else:
		email = sys.argv[1]
		api_key = sys.argv[2]
		in_dir = sys.argv[3]

		#Entrez.email = "arnab22.iitkgp@gmail.com"
		#Entrez.api_key = "51efc5e252e63fae1155a67c802bdd8e3e09"

		sp_all_tup = species_select(in_dir)
		Entrez.email = email
		Entrez.api_key = api_key

		species_iterate(sp_all_tup[0], 'id_list_major')
		species_iterate(sp_all_tup[1], 'id_list_large')
		species_iterate(sp_all_tup[2], 'id_list_interm')
		species_iterate(sp_all_tup[3], 'id_list_minor')

if __name__ == "__main__":
	main()
