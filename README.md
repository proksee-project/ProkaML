# ProkaML  
ProkaML is a machine learning based framework for prokaryotic sequence assembly quality control (QC) evaluation. The key steps in ProkaML are: building a database of prokaroyotic species with different genome assembly metrics, preprocessing and normalizing assembly attributes, and using random forest classification on the normalized genomic attributes to generate a model for assembly QC evaluation.  
This document provides an overview of the programs and steps used to accomplish the aforementioned steps. Additional details are provided in `WORKFLOW.md`.    

## Database building  
The scripts for generating proksee database comprising genomic attributes of NCBI contig assemblies are in the directories `database_build` and `add_genomic_attributes`. The database can be built using `Snakemake` workflow using a single core or in a cluster environment (for upto 10 jobs in parallel)

Usage: 
```
cd database_build   
snakemake --cores 1 --config email="dummy@email.com" --config api_key="dummy_api_key_01234"   
```  

where dummy email and API key entries should be replaced with user specific actual values. 
For running snakemake in a cluster environment, replace the second line of the former command with:  

```  
snakemake -j 10 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -o {cluster.output} -e {cluster.error}" --config email="dummy@email.com" --config api_key="dummy_api_key_01234"
```  

where `cluster.json` is a file with parameter specifications for sbatch command. `{cluster.partition}` is currently specified as `dummy_partition` and must be replaced with the actual partition name within the cluster. 

### Step-wise python scripts
Step 1: Run `get_assembly_UIDs.py` . Usage: 
```
python get_assembly_UIDs.py email api_key
```

Runs API queries on NCBI to obtain assembly unique identifiers (UIDs), which are written in batches of 10,000 to different files bearing format `entrez_id_list/Assembly_UID_chunk*.txt` where `*` is an integer  

Step 2: Run `get_entrez_metadata.py` . Usage: 
```
python get_entrez_metadata.py email api_key entrez_id_list/Assembly_UID_chunk{i}.txt
```

For a given file containing assembly UIDs, this script obtains different genomic assembly attributes from NCBI using biopython's esummary function. Assembly attributes are written in tab separated columns for every row UID to the directory `entrez_species_metadata`. Output files will have names in the format of `entrez_species_metadata/Assembly_UID_chunk{i}_metadata.txt`.  

Step 3: Run `log_entrez_metadata.py`. Usage: 
```
python log_entrez_metadata.py
```  

This examines the metadata files from Step 2 and logs the number of UIDs that were successfully queried to obtain metadata and number of UIDs that failed as well. For every file, a temporary log file of the format `log_entrez_metadata_chunk{i}.txt` is generated and upon generating all log files, aggregates are calculated and the log files are deleted.  

Step 4: Run `get_species_counts.py`. Usage: 
```
python get_species_counts.py email api_key
```

The metadata generated in Step 2 is organized in a species specific manner. For every species, a lowerbound threshold of 10 assemblies is applied for inclusion. Every species is also queried against the NCBI taxonomy database and only those corresponding to prokaryotes (taxonomy kingdom either Bacteria or Archaea) are included. Taxonomy functions written within `get_species_taxonomy.py` are imported to the current script. Output files are written in the format `species_reorganized_metadata/{j}_metadata.txt` where `{j}` corresponds to a species name joined by underscore (e.g. `species_reorganized_metadata/Salmonella_enterica_metadata.txt`)  

Step 5: Run `append_additional_attributes.py`. Usage: 
```
python append_additional_attributes.py email api_key species_reorganized_metadata/{j}_metadata.txt
```

This takes a species metadata file from step 4 and for every assembly row, downloads the complete fasta assembly, calculates the overall GC content and appends an additional column for every row. Files are written to the directory `additional_species_metadata` with names formatted to `{j}_added_attributes.txt`  

Step 6: Run `log_gc_content.py`. Usage: 
```
python log_gc_content.py
```  

This examines the metadata files with appended GC content from step 5 and logs the number of assemblies for which GC content was successfully calculated and also the numbers of assemblies that failed. For every file, a temporary log file of the format `{j}_log_gc_content.txt` is generated, where {j} corresponds to a species' name. Upon generating all log files, aggregates are calculated and the log files are deleted.  

Step 7: Run `concatenate_metadata.py`. Usage: 
```
python concatenate_metadata.py
```

This step concatenates all metadata files to generate four output files.  
Metadata for species with at least 1000 assemblies each are written to `major_species_metadata_added_attributes.txt`.  
Metadata for species having at least 100 assemblies each but not exceeding 1000 assemblies are written to `large_species_metadata_added_attributes.txt`.  
Metadata for species having at least 10 assemblies each but not exceeding 100 assemblies each are written to `intermediate_species_metadata_added_attributes.txt`.  
All metadata files are concatenated in a resulting file `well_represented_species_metadata_added_attributes.txt` containing assembly metadata for species with at least 10 assemblies. This file is used in subsequent steps of preprocessing, normalizing and generating machine learning models.
