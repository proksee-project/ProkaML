import constants as const
import glob
import os

log = open(const.LOG_FILE + const.FILE_EXTENSION, mode=const.APPEND_MODE)
log.write('\n#########################################################\n')
log.write('Getting NCBI metadata from Entrez esummary API query\n')
log.write('#########################################################\n')

total_entrez_metadata_success_UIDs = 0
total_entrez_metadata_failure_UIDs = 0

for chunk_log_file in glob.glob('log_entrez_metadata_chunk*.txt'):
	with open(chunk_log_file) as f:
		report = f.read().rstrip().split('\t')
		total_entrez_metadata_success_UIDs += int(report[1])
		total_entrez_metadata_failure_UIDs += int(report[2])
		if int(report[2]) > 0:
			message = const.UID_PREFIX + report[0] + const.FILE_EXTENSION + ' : ' + report[2] + " UID's metadata could not be obtained: \n"
			log.write(message)
			log.write(report[3] +'\n')

	os.remove(chunk_log_file)

success_message = 'Metadata obtained for a total number of ' + str(total_entrez_metadata_success_UIDs) + ' UIDs\n'
log.write(success_message)
failure_message = 'Metadata could not be obtained for a total number of ' + str(total_entrez_metadata_failure_UIDs) + ' UIDs\n'
log.write(failure_message)
