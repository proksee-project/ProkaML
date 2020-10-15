import os
import sys
from Bio import Entrez
import urllib.request
import urllib.error
from metadata_class import Metadata

#outfile.write('Organism Name\tOrganism Detail\tAssembly Name\tGenbank Accession\tRefseq Accession\tGenome Coverage\tSubmission Date\tLast Update Date\tRefseq Exclusion Reason\tContigN50\tAssembly Method\tSequencing Technology\n')

def idlist_parser(in_dir=None, out_dir=None, filename=None):
	with open(os.path.join(in_dir, filename)) as f:
		idlist = f.read().splitlines()
	
	output_filename = filename.split('idlist.txt')[0] + 'metadata.txt'
	output_file = open(os.path.join(out_dir, output_filename), 'w')
	
	return idlist, output_file

def main():
	input_dir = sys.argv[1]
	output_dir = sys.argv[2]
	filename = sys.argv[3]

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	idlist, output_file = idlist_parser(input_dir, output_dir, filename)

	Entrez.email = "arnab22.iitkgp@gmail.com"
	Entrez.api_key = "51efc5e252e63fae1155a67c802bdd8e3e09"

	metadata = Metadata(idlist, output_file)
	metadata.meta_integrate()

if __name__ == "__main__":
	main()




