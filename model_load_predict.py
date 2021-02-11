import pandas as pd
import os
from collections import defaultdict
import numpy as np

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
def incoming_dataframe_parser():
	input_df = pd.read_table('test_input_ml.txt')

	return input_df

#Currently want to check if species in input file are mapped correctly to their normalized metrics from dictionary
#Need to normalize metrics (calculate log, perform subtraction) for the species in separate functions to be written 
def main():
	sp_dicn = median_log_reader()

	input_df = incoming_dataframe_parser()
	sp_list = input_df['Organism Name'].to_list()
	for i in range(0, len(sp_list)):
		print(sp_list[i],'\t',sp_dicn[sp_list[i]])


if __name__ == '__main__':
	main()