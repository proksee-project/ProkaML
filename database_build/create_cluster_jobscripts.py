import os
import glob
from collections import defaultdict
import argparse


my_parser = argparse.ArgumentParser(usage='python %(prog)s [-h] email api_key output_dir partition_name',
                                    description='Creates scripts for cluster submission')
my_parser.add_argument('email',
                        type=str,
                        help='user email address')
my_parser.add_argument('api_key',
                        type=str,
                        help='NCBI user API key')
my_parser.add_argument('output_dir',
                        type=str,
                        help='output directory for cluster jobs')
my_parser.add_argument('partition_name',
                        type=str,
                        help='partition name for cluster jobs')						
args = my_parser.parse_args()

email = args.email
api_key = args.api_key
output_dir = args.output_dir
partition = args.partition_name

num_files = defaultdict(int)
category = ['major', 'large', 'interm']

for i in range(0, len(category)):
	num_files[category[i]] = len(glob.glob('id_list_' + category[i] + '/*_idlist.txt'))

line1 = '#!/bin/bash' + '\n'
line2 = '#SBATCH --job-name=biopython' + '\n'
line3 = '#SBATCH --output=' + output_dir + 'job_%J_out.txt' + '\n'
line4 = '#SBATCH --error=' + output_dir + 'job_%J_err.txt' + '\n'
line5 = '#SBATCH --partition=' + partition + '\n\n'

for key, value in num_files.items():
	jobscript_filename = open('jobscript_' + key + '.sh', 'w')
	command_filename = open('command.sh', 'a')

	jobscript_filename.writelines([line1, line2, line3, line4, line5])
	line6 = 'python metadata_print_fileindex.py ' + email + ' ' + api_key + ' ' + \
		'id_list_' + key + ' ' + key + '_species_metadata $SLURM_ARRAY_TASK_ID' + '\n'
	jobscript_filename.write(line6)

	if value > 1000:
		max_value = 1000 
	else:
		max_value = value

	if key == 'major':
		array_throttle = str(9)
	else:
		array_throttle = str(1)

	command_line = 'sbatch --array [0-' + str(max_value) + ']%' + array_throttle + ' jobscript_' + key + '.sh' + '\n'
	command_filename.write(command_line)

