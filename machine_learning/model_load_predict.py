import pandas as pd
import os
from collections import defaultdict
import numpy as np
import joblib

#load species median log metrics as dictionary
def median_log_reader():
	reader_fh = open('species_median_log_metrics.txt', 'r')
	sp_metrics = defaultdict(list)
	for line in reader_fh:
		row = line.rstrip().split('\t')
		sp_info = []
		
		'''
		species specific median logn50, logcontigcount, logl50, logtotlen, logcoverage 
		written in order
		'''
		for element in range(1, len(row)):
			try:
				sp_info.append(float(row[element]))
			except ValueError:
				pass
		sp_metrics[row[0]] = sp_info

	return sp_metrics

#load tsv input data to be evaluated
#currently customized as pandas dataframe with columns for species, n50, contig count, l50, totlen and coverage 
#will need to input these from quast output 
def incoming_dataframe_parser(species_median_dicn):
	species = 'Neisseria meningitidis'
	coverage = 100
	n50 = 5373
	contig_count = 1128
	l50 = 117
	totlen = 2393084

	try:
		input_logcoverage = round(np.log10(coverage),3)
		normalized_coverage = input_logcoverage - species_median_dicn[species][4]

		input_logn50 = round(np.log10(n50),3)
		normalized_n50 = input_logn50 - species_median_dicn[species][0]

		input_logcontigcount = round(np.log10(contig_count),3)
		normalized_contigcount = input_logcontigcount - species_median_dicn[species][1]

		input_logl50 = round(np.log10(l50),3)
		normalized_l50 = input_logl50 - species_median_dicn[species][2]

		input_logtotlen = round(np.log10(totlen),3)
		normalized_totlen = input_logtotlen - species_median_dicn[species][3]

		X_test = [normalized_n50, normalized_contigcount, normalized_l50, normalized_totlen, normalized_coverage]
		X_test_input = np.reshape(X_test, (1, -1))

		X_test_minus_coverage = [normalized_n50, normalized_contigcount, normalized_l50, normalized_totlen]
		X_test_input_minus_coverage = np.reshape(X_test_minus_coverage, (1, -1))
	
	except IndexError:
		#Species dictionary for its median metrics does not exist
		pass

	return (X_test_input, X_test_input_minus_coverage)

def random_forest_prediction(X_test, X_test_minus_coverage):
	loaded_rf1 = joblib.load('random_forest_n50_contigcount_l50_totlen_coverage.joblib')
	pred_nparr1 = loaded_rf1.predict_proba(X_test)
	pred_val1 = pred_nparr1[0,0]
	
	loaded_rf2 = joblib.load('random_forest_n50_contigcount_l50_totlen.joblib')
	pred_nparr2 = loaded_rf2.predict_proba(X_test_minus_coverage)
	pred_val2 = pred_nparr2[0,0]

	return pred_val1, pred_val2

#Currently want to check if species in input file are mapped correctly to their normalized metrics from dictionary
#Need to normalize metrics (calculate log, perform subtraction) for the species in separate functions to be written
#Need to load random forest model (random_forest_model.joblib) 
def main():
	sp_dicn = median_log_reader()
	X_test, X_test_minus_coverage = incoming_dataframe_parser(sp_dicn)
	print(X_test)
	pred_prob1, pred_prob2 = random_forest_prediction(X_test, X_test_minus_coverage)
	print(pred_prob1, pred_prob2)

if __name__ == '__main__':
	main()