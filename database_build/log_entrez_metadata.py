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

import constants as const
import glob
import os
import re


def main():
	log = open(const.LOG_FILE + const.FILE_EXTENSION, mode=const.APPEND_MODE)
	log.write('\n#########################################################\n')
	log.write('Getting NCBI metadata from Entrez esummary API queries\n')
	log.write('#########################################################\n')

	total_entrez_metadata_success_UIDs = 0
	total_entrez_metadata_failure_UIDs = 0
	failed_UIDs_pattern = re.compile(r'\[(.*)\]')
	list_failure_logs = []

	for chunk_log_file in glob.glob('log_entrez_metadata_chunk*.txt'):
		with open(chunk_log_file) as f:
			report = f.read().rstrip().split('\t')
			total_entrez_metadata_success_UIDs += int(report[1])
			failed_UIDs = failed_UIDs_pattern.search(report[2])
			if len(failed_UIDs.group(1)) > 0:
				num_failed_UIDs = len(failed_UIDs.group(1).split(','))
				total_entrez_metadata_failure_UIDs += num_failed_UIDs
				message = const.UID_OUTPUT_DIR + '/' + const.UID_PREFIX + report[0] + const.FILE_EXTENSION + \
					' : ' + str(num_failed_UIDs) + " UID's metadata could not be obtained: \n" + report[2] + '\n'
				list_failure_logs.append(message)

		os.remove(chunk_log_file)

	success_message = 'Metadata obtained for a total number of ' + str(total_entrez_metadata_success_UIDs) + ' UIDs\n'
	log.write(success_message)
	if total_entrez_metadata_failure_UIDs > 0:
		failure_message = 'Metadata could not be obtained for a total number of ' + str(total_entrez_metadata_failure_UIDs) + \
			' UIDs\n'
		log.write(failure_message)
		log.write('See corresponding file names and UIDs for which metadata could not be obtained\n')
		for failure_log in list_failure_logs:
			log.write(failure_log)

if __name__ == '__main__':
	main()
