#!/bin/bash
#SBATCH --job-name=biopython
#SBATCH --output=/home/CSCScience.ca/amandal/proksee-database/snakemake/database_build/job_%J_out.txt
#SBATCH --error=/home/CSCScience.ca/amandal/proksee-database/snakemake/database_build/job_%J_err.txt
#SBATCH --partition=NMLResearch

python metadata_print_fileindex.py arnab22.iitkgp@gmail.com  51efc5e252e63fae1155a67c802bdd8e3e09 id_list_major major_species_metadata $SLURM_ARRAY_TASK_ID
