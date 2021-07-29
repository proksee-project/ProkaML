#!/bin/bash
#SBATCH --job-name=biopython
#SBATCH --output=dummy_directoryjob_%J_out.txt
#SBATCH --error=dummy_directoryjob_%J_err.txt
#SBATCH --partition=NMLResearch

python metadata_print_fileindex.py dummy@email.com dummy_api_key_01234 id_list_large large_species_metadata $SLURM_ARRAY_TASK_ID
