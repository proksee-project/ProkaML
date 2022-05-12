import numpy as np
import pandas as pd

df = pd.read_csv('well_represented_species_metadata_top100.txt', sep='\t')

def split_dataframe(df, chunk_size = 10): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

df_chunks = split_dataframe(df)

for i in range(len(df_chunks)):
	chunk_filename = 'assembly_annotated_metadata_' + str(i+1) + '.txt'
	df_chunks[i].to_csv(chunk_filename, sep='\t', index=False)

