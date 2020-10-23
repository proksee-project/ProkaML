# proksee-database
Reference library for microorganisms with different genome assembly metrics and statistics

Step 1: Run `assemblydb_entrez_query.py` . This runs Entrez API queries to scan for contig assemblies on entire NCBI assembly database. Generates counts of species/organism names in a two column tab separated text file `species_counts.txt`  

Step 2: Run `idlist_retriever.py`. This examines `species_counts.txt` and generates list of assembly UIDs for all species. Since some species are better represented in assembly database than others, this script writes UIDs into four separate folders based on counts of assemblies for a particular species  
- id_list_major - Folder containing UIDs for species >= 1000 assemblies. Species with > 10,000 assemblies are separated into chunk files each containing maximum of 10,000 UIDs each    
- id_list_large - Folder containing UIDs for species < 1000 but >= 100 assemblies each   
- id_list_interm - Folder containing UIDs for species < 100 but >= 10 assemblies each
- id_list_minor - Folder containing UIDs for species < 10 assemblies each
