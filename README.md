# proksee-database
proksee database is a reference library for microorganisms with different genome assembly metrics and statistics.  

The scripts for generating proksee database comprising genomic attributes of NCBI contig assemblies are in the directory `database_build`. The scripts run biopython API queries on the NCBI database and therefore require the user to provide an API key corresponding to their account. If you dont have an NCBI API key, visit the NCBI login [page](https://www.ncbi.nlm.nih.gov/account/) and create an account using your email. Once you sign in, click on the top right corner on your email ID. This will redirect you to a new page where you can find your API key. You will need to use this API key for most of the scripts within `database_build` directory.

## proksee-database generation (Snakemake automated workflow)  
Usage: `snakemake --cores 1 --config email="dummy@email.com" --config api_key="dummy_api_key_01234"`.  
where dummy email and API key entries should be replaced with user specific actual values.  
The generation of **proksee-database** utilizes snakemake workflow to automate complex processes with a single command. It is advisable to run snakemake in a separate conda environment. Visit the snakemake [page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for instructions on installation and setting up the appropriate software environment.  
In summary, the snakemake file `Snakefile` chains together individual python scripts, with their respective input and output dependencies to generate **proksee-database**. Currently, the database build time is approximately 10 days, which can be reduced by implementing parallel workflows in a cluster environment. To run snakemake in a cluster environment, user must have access to a cluster and the following command can be run:

`snakemake -j 10 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -t {cluster.time}" --config email="dummy@email.com" --config api_key="dummy_api_key_01234"` 


For understanding the operations of individual python scripts, users are requested to refer to the section below.   

## proksee-database generation (step-wise scripts)
Step 1: Run `assemblydb_entrez_query.py` . Usage: `python assemblydb_entrez_query.py email api_key`. 
This runs Entrez API queries to scan for contig assemblies on the entire NCBI assembly database. Generates counts of species/organism names in a two column tab separated text file `species_counts_[Month]_[year].txt`  

Troubleshoot: If the program crashes with the following error message:
```  
species = docsum['DocumentSummarySet']['DocumentSummary'][j]['SpeciesName']
IndexError: list index out of range
```    
This is possibly due to some server issue in NCBI. The script can be re-run without any errors (hopefully).  
Expected output: Two column tab separated text file `species_counts_[Month]_[year].txt`. This may look like:  
```  	
Salmonella enterica	281497
Escherichia coli	85549
Campylobacter jejuni	38238
Listeria monocytogenes	34181
Campylobacter coli	15398
.........................
.........................  
```  

Step 2: Run `idlist_retriever.py`. Usage: `python assemblydb_entrez_query.py email api_key`. 
This examines `species_counts_[Month]_[year].txt` from step 1 and generates list of assembly UIDs for all species. Since some species are better represented (more assembly counts) in assembly database than others, this script writes UIDs into four separate directories based on counts of assemblies for a particular species  
- id_list_major - Directory containing UIDs for species >= 1000 assemblies each. Species with > 10,000 assemblies are separated into chunk files each containing maximum of 10,000 UIDs      
- id_list_large - Directory containing UIDs for species < 1000 but >= 100 assemblies each   
- id_list_interm - Directory containing UIDs for species < 100 but >= 10 assemblies each
- id_list_minor - Directory containing UIDs for species < 10 assemblies each

Step 3 (Recommended): Run `metadata_print_fileindex.py`. This takes in one of the input directories generated in Step 2, captures NCBI metadata corresponding to different genomic assemmbly metrics and writes the metadata for all files containing assembly UIDs to a desired output directory. Depending on the number of files in a directory, the script takes integer index corresponding to the file number as an argument. This facilitates parallel processing.  
e.g. if `id_list_interm` dirctory has 10 files, and the desired output directory is `interm_metadata`, the script `metadata_print_fileindex.py` can be used in the following way:  
- `python metadata_print_fileindex.py email api_key id_list_interm interm_metadata 0`  
- `python metadata_print_fileindex.py email api_key id_list_interm interm_metadata 1`  
- `python metadata_print_fileindex.py email api_key id_list_interm interm_metadata 2`  
- `.................................................................................`  
- `python metadata_print_fileindex.py email api_key id_list_interm interm_metadata 9`  

The parallel processing should be performed in a cluster or a computer with sufficient number of available cores. In theory, `metadata_print_fileindex.py` can be run in as many parallel instances depending on the number of files, but care must be taken so as not to run more than 10 parallel instances. This is because NCBI limits the number of API requests to 10 per second. Over-usage of NCBI API requests may block the user or even an IP of an institution. Should the need for parallel extensive API requests arise, user should contact eutilities@ncbi.nml.nih.gov.

Step 3 (Alternate): Run `metadata_print_filename.py`. This step is almost identical to the recommended step 3 with the only difference arising in providing filename as the last argument instead of life number index.  
e.g. if one of the files within directory `id_list_interm` is `Acholeplasma_laidlawii_chunk1_idlist.txt`, the usage of this script would be: 
- `python metadata_print_filename.py email api_key id_list_interm interm_metadata Acholeplasma_laidlawii_chunk1_idlist.txt`  

The output file for the example would be `Acholeplasma_laidlawii_chunk1_metadata.txt`. The alternate step is used for troubleshooting when individual files containing UIDs are not annotated with NCBI metadata for reasons beyond the user's control (NCBI server/connection issues).  

Step 3 (either recommended or alternate) writes a tab separated table `*chunk*_metadata.txt` with rows corresponding to every assembly within a `*chunk*_idlist.txt` file and columns corresponding to the following genomic attributes:  
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

The genomic attributes are obtained from the `get_genomic_metadata.py` script from the `AttributeMetadata` class. A file named `Acholeplasma_laidlawii_chunk1_idlist.txt` containing UIDs will generate `Acholeplasma_laidlawii_chunk1_metadata.txt` upon the running of Step 3.

Step 4: The annotated metadata files are concatenated. The naming of concatenated files is based on the grouping in Step 2. *E.g.* species with >= 1000 assembly records, with UIDs in `id_list_major` and with different metadata annotated fileparts (from Step 3) are concatenated as `major_species_metadata.txt`. Overall, there are four concatenated metadata files named as:  
- `major_species_metadata.txt` for species >= 1000 assemblies each  
- `large_species_metadata.txt` for species < 1000 but >= 100 assemblies each  
- `intermediate_species_metadata.txt` for species < 100 but >= 10 assemblies each  
- `minor_species_metadata.txt` for species < 10 assemblies each  
In order to conduct analytical strategies with well represented species, the `minor_species_metadata.txt` file containing metadata for species with poor representation in the NCBI assembly database is not used. The other annotated metadata files: `major_species_metadata.txt`, `large_species_metadata.txt` and `intermediate_species_metadata.txt` are concatenated together as `well_represented_species_metadata.txt` with a one-line column header file `metadata_header.txt` appended at the beginning. The resulting file forms the starting point for subsequent analyses.  
