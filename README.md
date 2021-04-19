# proksee-database
Reference library for microorganisms with different genome assembly metrics and statistics

## proksee-database generation 
The scripts for generating proksee database comprising genomic attributes of NCBI contig assemblies are in the directory `database_build`. The scripts run biopython API queries on the NCBI database and therefore require the user to provide an API key corresponding to their account. If you dont have an NCBI api key, create an account at https://www.ncbi.nlm.nih.gov/account/ using your email. Once you sign in, click on the top right corner on your email id. This will redirect you to a new page where you can find your api key. You will need to use this api key for most of the scripts within `database_build` directory.  

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

The output file for the example would be `Acholeplasma_laidlawii_chunk1_metadata.txt`. The alternate step is used for troubleshooting when individual files containing UIDs are not annotated with NCBI metadata for reasons beyond the user's control (NCBI server/connection issues)  

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

The genomic attributes are obtained from `GetMetadata` class in `get_genomic_attributes.py`. A file `Acholeplasma_laidlawii_chunk1_idlist.txt` containing UIDs will generate `Acholeplasma_laidlawii_chunk1_metadata.txt` upon running of Step 3.

Step 4: The annotated metadata files are concatenated. The naming of concatenated files are based on the grouping in step 2. e.g. species with > 1000 assembly records, with UIDs in `id_list_major` and with different metadata annotated fileparts (from step 3) are concatenated as `major_species_metadata.txt`. Overall, there are four concatenated metadata files named as:  
- `major_species_metadata.txt` for species >= 1000 assemblies each  
- `large_species_metadata.txt` for species < 1000 but >= 100 assemblies each  
- `intermediate_species_metadata.txt` for species < 100 but >= 10 assemblies each  
- `minor_species_metadata.txt` for species < 10 assemblies each  
In order to conduct analytical strategies with well represented species, the `minor_species_metadata.txt` file containing metadata for species with poor representation in the NCBI assembly database is not used. The other annotated metadata files : `major_species_metadata.txt`, `large_species_metadata.txt` and `intermediate_species_metadata.txt` are concatenated together as `well_represented_species_metadata.txt` with a one-line column header file `metadata_header.txt` appended at the beginning. The resulting file forms the starting point for subsequent analyses.  

## adding additional genomic attributes  
The scripts for adding custom genomic attributes are in the directory `add_genomic_attributes`.  
Usage: `python add_genomic_attributes.py EMAIL API_KEY FILE_NAME OUTPUT_DIRECTORY`  
This program takes the metadata file generated from the final step of **proksee-database generation**, computes genomic attributes of interest for every assembly and writes the genomic attributes as additional columns to an output spreadsheet. For example, the file `well_represented_species_metadata.txt` will be processed and re-written as `well_represented_species_metadata_added_attributes.txt`. Other inputs for the script `add_genomic_attributes.py` include email address, NCBI api key and a user-desired output directory for downloading intermediate files. Currently, the only additional calculated attribute is overall GC content, which is the fraction of G or C bases of all nucleotide bases of an assembly. GC content is calculated from `gc_content.py`. Other genomic attributes of interest are in development phase and may be added later. 

## preprocessing and normalizing  
`clean_metadata.py` organizes assembly methods and sequencing technology using regular expressions. Subsequently rows for long read data are identified from assembly methods and sequencing technology and excluded.   
`preprocess_metadata.py` performs species specific grouped normalizing of metadata, firstly by calculating species specific medians (or medians of logarithms) of genomic attributes of interest and then normalizing the genomic attributes with respect to their species specific median computed values.  

## building and exporting machine learning model
`cmd_machine_learning.py` subsets the dataframe with curated labels (RefSeq included/excluded) and fits random forest classifier models on curated data. Models are evaluated over a 10 fold cross validation and the best performing model is output as a joblib object for subsequent analyses
