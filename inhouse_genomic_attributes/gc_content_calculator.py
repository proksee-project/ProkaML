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
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
from Bio import Entrez 
import urllib.request
import urllib.error
import sys
import gzip
import re
import time
import zlib

"""
Given a spreadsheet of NCBI assembly accession IDs (or Genbank accession IDs), this program 
downloads the contig assemblies from NCBI ftp server, calculates the overall GC content as 
a fraction of the total length, and writes the GC content as an additional column 
"""

#Use pandas to generate a list of Genbank accession IDs from an input spreadsheet
def metadata_file_read(file):
	df = pd.read_csv(os.path.join(file), sep='\t')
	gb_acc = df['Genbank Accession'].to_list()

	return df, gb_acc

#Use biopython Entrez library to obtain ftp downlo
def assem_dnld(gb_id, out_dir):
	try:
		#Generate ftp download link of contig assembly 
		handle = Entrez.esearch(db="assembly", term=gb_id)
		record = Entrez.read(handle)
		handle2 = Entrez.esummary(db="assembly", id=record['IdList'])
		record2 = Entrez.read(handle2, validate=False)
		gb_dir = record2['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
		seq_file = os.path.basename(gb_dir) + '_genomic.fna.gz'
		seq_file_link = os.path.join(gb_dir, seq_file)
		
		#provide path in to download the assembly file
		full_file_path = os.path.join(out_dir, seq_file)

		#Avoid duplicate downloading if program interrupts
		if download_needed(seq_file_link, full_file_path):
			urllib.request.urlretrieve(seq_file_link, full_file_path)
	
	except Exception:
		full_file_path = 'NA'

	return full_file_path

#Check for file size to determine if assembly file already exists
def download_needed(url_link, full_file_path):
	url_request = urllib.request.Request(url_link, method='HEAD')
	url_open = urllib.request.urlopen(url_request, timeout=100)
	file_size = int(url_open.headers['Content-Length'])
	if (not os.path.exists(full_file_path) or \
			file_size != os.path.getsize(full_file_path)):
		return True
	else:
		return False

#Calculate gc content as a fraction of the total assembly length 
def calculate_gc_func(full_file_path):
	if full_file_path == 'NA':
		overall_gc_content = 'NA'
	else:
		gc_content = 0
		full_length = 0
		try:
			with gzip.open(full_file_path, mode='rt') as open_file:
				for line in open_file:
					if not line.startswith('>'):
						full_length += len(line.rstrip('\n'))
						nt_uppercase = line.rstrip('\n').upper()
						for gc_pattern in 'GC':
							gc_content += nt_uppercase.count(gc_pattern)
			overall_gc_content = round(gc_content/full_length, 3)
		except Exception:
			overall_gc_content = 'NA'
	return overall_gc_content

def main():
	filename = sys.argv[1]
	out_dir = sys.argv[2]
	df, gb_acc_list = metadata_file_read(filename)
	
	#Use your real email and corresponding NCBI api key
	Entrez.email = "youremail@email.com"
	Entrez.api_key = "yourapikey"
	
	#Iterate through list of Genbank IDs to download assembly file and calculate GC content
	gc_content_list = []
	for i in range(0, len(gb_acc_list)):
		file_path = assem_dnld(gb_acc_list[i], out_dir)
		print('{} assembly downloaded.'.format(gb_acc_list[i]), end='')
		
		#append the gc content to a list
		gc_content_list.append(calculate_gc_func(file_path))
		print('gc content calculated.', sep='')
		
		"""
		Remove the downloaded file. This block can be commented out if assemblies do not take
		up much of hard drive space
		"""
		try:
			os.remove(file_path)
			print('{} assembly file removed'.format(gb_acc_list[i]))
		except FileNotFoundError:
			pass

			#Write the list of GC content as a column to output file
	df['GCcontent'] = gc_content_list
	outfile = filename.split('.')[0] + '_gccalculated.txt'
	outpath = os.path.join(out_dir, outfile)
	df.to_csv(outpath, sep='\t', mode='w', index=False)

if __name__ == "__main__":
	main()
