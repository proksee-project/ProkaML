# proksee-database
Reference library for microorganisms with different genome assembly metrics and statistics

Step 1: Run `assemblydb_entrez_query.py` . This runs Entrez API queries to scan for contig assemblies on entire NCBI assembly database. Generates counts of species/organism names in a two column tab separated text file `species_counts.txt`  

Step 2: Within the `metadata` folder, run `idlist_retriever.py`. This examines `species_counts.txt` and generates list of assembly UIDs for all species. Since some species are better represented in assembly database than others, this script writes UIDs into four separate folders based on counts of assemblies for a particular species  
- id_list_major - Folder containing UIDs for species >= 1000 assemblies each. Species with > 10,000 assemblies are separated into chunk files each containing maximum of 10,000 UIDs      
- id_list_large - Folder containing UIDs for species < 1000 but >= 100 assemblies each   
- id_list_interm - Folder containing UIDs for species < 100 but >= 10 assemblies each
- id_list_minor - Folder containing UIDs for species < 10 assemblies each

Step 3: Within the `metadata` folder, run `metadata_extract_fileindex.py`. This takes in one of the input directories generated in Step 2, captures NCBI metadata corresponding to different genomic assemmbly metrics and writes the metadata for all files containing assembly UIDs to a desired output directory. Depending on the number of files in a directory, the script takes integer index as an argument which facilitates parallel processing.  
e.g. if id_list_interm dirctory has 10 files, and the desired output directory is interm_metadata, the script `metadata_extract_fileindex.py` can be used in the following way:
metadata_extract_fileindex.py id_list_interm interm_metadata 0
metadata_extract_fileindex.py id_list_interm interm_metadata 1
metadata_extract_fileindex.py id_list_interm interm_metadata 2
..............................................................
metadata_extract_fileindex.py id_list_interm interm_metadata 9

The parallel processing should be performed in a cluster or a computer with sufficient number of available cores. In theory, `metadata_extract_fileindex.py` can be run in as many parallel instances depending on the number of files, but care must be taken so as not to run more than 10 parallel instances. This is because NCBI limits the number of API requests to 10 per second. Overusage may block the user or even an IP of an institution. Should the need for parallel extensive API requests arise, user should contact eutilities@ncbi.nml.nih.gov.

Step 4: 
