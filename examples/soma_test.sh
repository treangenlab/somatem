#! /bin/bash

#SBATCH --job-name=somatem
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --time=23:30:00
#SBATCH --mail-user=mail@mail.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --threads-per-core=1
#SBATCH --mem=250GB
#SBATCH --export=ALL

source /path/to/miniforge3/etc/profile.d/conda.sh
conda activate somatemtest

nextflow run /path/to/Documents/SOMAteM/subworkflows/somatem_mags.nf \
  --input_dir   /path/to/Documents/SOMAteM/examples/data/input4mags \
  --output_dir  /path/to/Documents/SOMAteM/examples/data/mag_output \
  --threads     96 \
  --flye_mode nano-hq \
  --semibin_environment human_gut \
  --checkm2_db  /path/to/checkm2/uniref100.KO.1.dmnd \
  --bakta_db    /path/to/bakta/db \
  -c /path/to/Documents/SOMAteM/confs/somatem_mags.config
