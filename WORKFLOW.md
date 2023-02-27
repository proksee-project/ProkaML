# ProkaML workflow  

## Database building  
The scripts in `database_build` and `add_genomic_attributes` directories run biopython API queries on the NCBI database and therefore require the user to provide an API key corresponding to their account. If you dont have an NCBI API key, visit the NCBI login [page](https://www.ncbi.nlm.nih.gov/account/) and create an account using your email. Once you sign in, click on the top right corner on your email ID. This will redirect you to a new page where you can find your API key. You will need to use this API key for most of the scripts within `database_build` directory.   
We utilize snakemake workflow to automate complex processes with a single command. It is advisable to run snakemake in a separate conda environment. Visit the snakemake [page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for instructions on installation and setting up the appropriate software environment.  
In summary, the snakemake file `Snakefile` chains together individual python scripts, with their respective input and output dependencies to generate **ProkaML database**. Currently, the database build time is approximately 20 days, which can be reduced by implementing parallel jobs in a cluster environment.  

Troubleshooting:  
Snakemake workflow in cluster requires a user to maintain active connection. If this is not a possibility, `sbatch` must be appended at the beginning of the command:  

`sbatch snakemake --use-conda -j 10 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -t {cluster.time}-o {cluster.output} -e {cluster.error}" --config email=="dummy@email.com" api_key="dummy_api_key_01234"`  

`sbatch` at the beginning ensures that the main snakemake workflow is always running, while `sbatch` at the middle lets snakemake parallelize the workflow. 

If snakemake cluster jobs have to be terminated or are killed for some reason, the following snakemake commands must be run before re-running snakemake instance:  
```
snakemake --cleanup-metadata <filenames>
snakemake --unlock
```

## Step-wise python scripts  
`get_species_assembly.py` runs Entrez API queries to scan for contig assemblies on the entire NCBI assembly database and generates counts of species/organism names in a two column tab separated text file `species_counts_[Month]_[year].txt`  

Troubleshooting:  
If the program crashes with the following error message:
```  
species = docsum['DocumentSummarySet']['DocumentSummary'][j]['SpeciesName']
IndexError: list index out of range
```    
This is possibly due to some server issue in NCBI and the script can be re-run. The correct output file `species_assemblycounts_[Month]_[year].txt` should look similar to:  
```  	
Salmonella enterica	281497
Escherichia coli	85549
Campylobacter jejuni	38238
Listeria monocytogenes	34181
Campylobacter coli	15398
.........................
.........................  
```  

`retrieve_idlist.py`. When number of assembly UIDs for a species exceeds 10,000, UIDs are written in species' specific separate chunk files, with each chunk containing not more than 10,000 UIDs.  

`get_entrez_metadata.py`. Imports `EntrezMetadata` class from `entrez_metadata.py`. The following genomic attributes are obtained for every assembly UID using Entrez esummary function:  
- Species name
- Species strain/isolate
- Assembly ID
- Genbank ID
- Refseq ID (or NA)
- Genome coverage
- Assembly submission date
- Assembly last update date
- Refseq Exclusion Reason (if Refseq ID is NA)
- N50
- Number of contigs
- L50
- Assembly length
- Assembly method
- Sequencing technology  

`append_additional_attributes.py`. Imports `GCContentCalculate` class from `gc_content.py` which is present in directory `add_genomic_attributes`. Full fasta assemblies are downloaded and deleted (post calculation of GC content) in an iterative manner and therefore, should not throttle hard drive storage.  

