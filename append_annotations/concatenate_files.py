import pandas as pd
import glob

df_list = []
for file in glob.glob('assembly_annotated_metadata_*.txt'):
	df = pd.read_csv(file, sep='\t')
	df_list.append(df)

concatenated_dataframe = pd.concat(df_list, axis=0)
concatenated_dataframe.to_csv('concat.txt', sep='\t', index=False)

