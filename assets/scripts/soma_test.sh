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

nextflow run /path/to/Somatem/subworkflows/somatem_mags_nf \
  --input_dir   /path/to/Somatem/examples/data/input4mags \
  --output_dir  /path/to/Somatem/examples/data/mag_output3 \
  --threads     96 \
  --flye_mode nano-hq \
  --semibin_environment human_gut \
  --completeness_threshold 50 \
  --checkm2_db  /path/to/checkm2/uniref100.KO.1.dmnd \
  --bakta_db    /path/to/bakta/db \
  --singlem_metapackage /path/to/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb \
  -c /path/to/SOMAteM/conf/simple_somatem_mags.config \
