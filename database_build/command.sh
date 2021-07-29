sbatch --array [0-98]%9 jobscript_major.sh
sbatch --array [0-215]%1 jobscript_large.sh
sbatch --array [0-1000]%1 jobscript_interm.sh
