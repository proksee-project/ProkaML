import constants as const
import glob
import os

log = open(const.LOG_FILE + const.FILE_EXTENSION, mode=const.APPEND_MODE)
log.write('\n#########################################################\n')
log.write('Getting NCBI metadata from Entrez esummary API query\n')
log.write('#########################################################\n')

total_gc_content_success= 0
total_gc_content_failures = 0

for sublog_file_list in glob.glob('*_log_gc_content.txt'):
	sublog_file = open(sublog_file_list)
	report = sublog_file.read().rstrip().split('\t')
	total_gc_content_success = int(report[0])
	total_gc_content_failures = int(report[1])

	os.remove(sublog_file)

success_message = 'GC content obtained for a total number of ' + str(total_gc_content_success) + ' assemblies\n'
log.write(success_message)
failure_message = 'GC content could not be obtained for a total number of ' + str(total_gc_content_failures) + ' assemblies\n'
log.write(failure_message)
