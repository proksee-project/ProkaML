'''
Copyright:

University of Manitoba & National Microbiology Laboratory, Canada, 2021

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

import pandas as pd
import subprocess
import re
import argparse
import os
import shutil
import glob

def main():
	my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] input_file',
                                    description='Maps genbank assemdly ID to assembly file')
	my_parser.add_argument('input_file',
                        type=str,
                        help='file containing assembly metadata')

	args = my_parser.parse_args()
	input_file = args.input_file
	output_file = open(input_file.split('.')[0] + '_buscostats.' + input_file.split('.')[1], 'w')

	line_count = 0
	with open(input_file, 'r') as assembly_lineage:
		for line in assembly_lineage:
			line_count += 1
			if line_count == 1:
				output_file.write(line.rstrip() + '\tCompleteness\tSingle BUSCOs\tDuplicate BUSCOs\n')

			else:
				row = line.rstrip().split('\t')
				genbank_id = row[3]
				assembly_filepath_zipped = 'additional_species_metadata/' + row[17]
				assembly_filepath_unzipped = assembly_filepath_zipped[:len(assembly_filepath_zipped) - 3]
				busco_lineage = row[18]

				failure = re.search(r'taxonomy not found', row[19])
				if failure is not None:
					output_file.write('\t'.join(row),'\tNA\tNA\tNA\n')
					print('\t'.join(row),'\tNA\tNA\tNA')

				else:
					busco_cmd = 'busco -i ' + assembly_filepath_unzipped + ' -o ' + genbank_id + ' -m genome -l ' + busco_lineage + ' -c 4 -f'
					print(busco_cmd)
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

					shutil.rmtree(genbank_id)
					shutil.rmtree('busco_downloads')
		
		busco_log_files = glob.glob('busco_*')
		for file in busco_log_files:
			os.remove(file)


if __name__ == '__main__':
	main()
