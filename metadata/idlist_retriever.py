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

def species_select(in_dir=None):
	query_sp_list = []
	infile = open(os.path.join(in_dir, 'species_counts.txt'),'r')
	for line in infile:
		sp = line.rstrip().split('\t')
		if int(sp[1]) > 1000:
			sp_tup = sp[0], sp[1]
			query_sp_list.append(sp_tup)

	return query_sp_list

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
	species_name_list = species_tuple[0].split()
	for i in range(1, chunks + 1):
		outfile = "_".join(species_name_list) + '_chunk' + str(i) + '_idlist.txt'
		species_outfile_list.append(outfile)

	return species_outfile_list

def species_iterate(species_list=None, out_dir=None):
	for sp_tup in species_list:
		idlist = id_list(sp_tup)
		chunks = chunk_counter(idlist)
		species_outfile_list = species_outfile(sp_tup, chunks)
		for k in range(0, chunks):
			idlist_lowerlim = k*10000
			idlist_upperlim = (k+1)*10000 + 1
			idlist_batch = []
			for m in range(idlist_lowerlim, min(idlist_upperlim, len(idlist))):
					idlist_batch.append(idlist[m])
			output_handle = os.path.join(out_dir, species_outfile_list[k])
			output_file = open(output_handle, 'w')
			idlist_string = '\n'.join(idlist_batch)
			output_file.write(idlist_string)

def main():
	if len(sys.argv) != 5:
		sys.exit('''
		Command usage: python idlist_retriever.py EMAIL NCBI_API_KEY INPUT_DIRECTORY OUTPUT_DIRECTORY. 
		Need to pass 3 arguments corresponding to your email, your ncbi api key and the input directory
		containing species_counts.txt
		''')
	
	else:
		email = sys.argv[1]
		api_key = sys.argv[2]
		in_dir = sys.argv[3]
		out_dir = sys.argv[4]

		if not os.path.exists(out_dir):
			os.mkdir(out_dir)
	
		#Entrez.email = "arnab22.iitkgp@gmail.com"
		#Entrez.api_key = "9b946150b1df235cf0e7c288ca4588a97c08"

		species_list = species_select(in_dir)
		Entrez.email = email
		Entrez.api_key = api_key

		species_iterate(species_list, out_dir)

if __name__ == "__main__":
	main()
