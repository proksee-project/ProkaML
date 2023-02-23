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
	"""
	Compiles reports of assembly UIDs for which metadata was sucessfully or unsuccessfully downloaded
	"""

	log = open(const.LogFiles.MAIN_LOG + const.FileFormat.TEXT , mode=const.FileFormat.APPEND_MODE)
	log.write('\n#########################################################\n')
	log.write('Downloading assemblies, appending GC content to the metadata\n')
	log.write('#########################################################\n')

	total_gc_success_ids = 0
	total_gc_fail_ids = 0
	failed_ids_pattern = re.compile(r'\[(.*)\]')
	list_failure_logs = []

	for chunk_log_file in glob.glob(os.path.join(const.FileDirectories.DATABASE_PATH, const.FileDirectories.SPECIES_GC_LOG_DIR, 
		'*' + const.LogFiles.SUB_LOG_GC + const.FileFormat.TEXT)):
		with open(chunk_log_file) as f:
			report = f.read().rstrip().split(const.FileFormat.SEPARATOR)
			total_gc_success_ids += int(report[0])
			failed_UIDs = failed_UIDs_pattern.search(report[2])
			if len(failed_UIDs.group(1)) > 0:
				num_failed_UIDs = len(failed_UIDs.group(1).split(','))
				total_entrez_metadata_failure_UIDs += num_failed_UIDs
				message = const.FileDirectories.UID_OUTPUT_DIR + '/' + const.Assembly.UID_PREFIX + report[0] + \
					const.FileFormat.FILE_EXTENSION + ' : ' + str(num_failed_UIDs) + " UID's metadata could not be obtained: \n" \
					+ report[2] + '\n'
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
	log.close()


if __name__ == '__main__':
	main()
