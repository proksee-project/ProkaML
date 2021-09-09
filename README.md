# ProkaML  
ProkaML is a machine learning based framework for prokaroyotic sequence assembly quality control (QC) evaluation. The key steps in ProkaML are: building a database of prokaroyotic species with different genome assembly metrics, preprocessing and normalizing assembly attributes, and using random forest classification on the normalized genomic attributes to generate a model for assembly QC evaluation.  
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
Step 1: Run `get_species_assemblies.py` . Usage: 
```
python get_species_assemblies.py email api_key
```

Generates counts of species/organism names in `species_assemblycounts_[Month]_[year].txt`  

Step 2: Run `get_taxonomical_kingdom.py` . Usage: 
```
python get_taxonomical_kingdom.py email api_key
```

This examines `species_assemblycounts_[Month]_[year].txt` from step 1 and queries every species to obtain its kingdom from the NCBI taxonomy database. Taxonomical kingdom names are appended to generate the file `species_assemblycounts_[Month]_[year]_taxonomy.txt`. 

Step 3: Run `retrieve_idlist.py`. Usage: 
```
python retrieve_idlist.py email api_key species_assemblycounts_[Month]_[year]_taxonomy.txt
```  

This examines `species_assemblycounts_[Month]_[year]_taxonomy.txt` Only species belonging to kingdoms Bacteria or Archaea are retained and assembly UIDs for all species are written in species specific separate files to the directory `entrez_id_list`.  

Step 4: Run `get_entrez_metadata.py`. Usage: 
```
python get_entrez_metadata.py email api_key entrez_id_list/{i}_idlist.txt
```

This takes assembly UIDs for a given file for a given species (denoted by `{i}`. As an example, `{i}` could be `Salmonella_enterica_chunk4`) and obtains different genomic assembly attributes for every UID from NCBI. Assembly attributes are written in tab separated columns for every row UID to the directory `entrez_species_metadata`. Output files will have names in the format of `{i}_metadata.txt`.  

Step 5: Run `append_additional_attributes.py`. Usage: 
```
python append_additional_attributes.py email api_key entrez_species_metadata/{i}_metadata.txt
```

This takes a metadata file from step 4 and for every assembly row, downloads the complete fasta assembly, calculates the overall GC content and appends an additional column for every row. Files are written to the directory `additional_species_metadat`a with names formatted to `{i}_metadata_added_attributes.txt`

Step 6: Run `concatenate_metadata.py`. Usage: `python concatenate_metadata.py`.  
This step concatenates all metadata files to generate four output files.  
Metadata for species with at least 1000 assemblies each are written to `major_species_metadata_added_attributes.txt`.  
Metadata for species having at least 100 assemblies each but not exceeding 1000 assemblies are written to `large_species_metadata_added_attributes.txt`.  
Metadata for species having at least 10 assemblies each but not exceeding 100 assemblies each are written to `intermediate_species_metadata_added_attributes.txt`.  
All metadata files are concatenated in a resulting file `well_represented_species_metadata_added_attributes.txt` containing assembly metadata for species with at least 10 assemblies. This file is used in subsequent preprocessing, normalizing and generating machine learning models.
