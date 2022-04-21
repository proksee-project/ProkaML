import pandas as pd
import subprocess
import re
import sys
filename = sys.argv[1]
import os
output_file = open(filename + '_buscostats.txt', 'w')
import shutil

with open(filename, 'r') as assembly_lineage:
	#next(assembly_lineage)
	for line in assembly_lineage:
		row = line.rstrip().split('\t')
		genbank_id = row[2]
		assembly_filepath = 'predownloaded_assemblies/' + row[5]
		busco_lineage = row[6]

		failure = re.search(r'taxonomy not found', row[7])
		if failure is not None:
			output_file.write('\t'.join(row),'\tNA\tNA\tNA\n')
			print('\t'.join(row),'\tNA\tNA\tNA')

		else:
			busco_cmd = 'busco -i ' + assembly_filepath + ' -o ' + genbank_id + ' -m genome -l ' + busco_lineage + ' -c 4 -f'
			with open('busco_' + genbank_id + '_out.txt', 'w+') as busco_out:
				subprocess.call(busco_cmd, shell=True, stdout=busco_out)
				busco_out.seek(0)
				for busco_line in busco_out.read().split('\n'):
					busco_stats = re.search(r'^.*?\|C:(.+?)%\[S:(.+?)%,D:(.+?)%\]', busco_line)
					if busco_stats is not None:
						completeness = busco_stats.group(1)
						single_buscos = busco_stats.group(2)
						duplicate_buscos = busco_stats.group(3)

						output_file.write('\t'.join(row) +'\t' + completeness + '\t' + single_buscos + '\t' + duplicate_buscos + '\n')
						print('\t'.join(row), '\t', completeness, '\t', single_buscos, '\t', duplicate_buscos)

		#shutil.rmtree(genbank_id)
		#busco_stdout_file = 'busco_' + genbank_id + '_out.txt'
		#os.remove(busco_stdout_file)