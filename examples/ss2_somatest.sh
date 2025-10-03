#! /bin/bash

#SBATCH --job-name=msomatem
#SBATCH --partition=debug
#SBATCH --account=commons
#SBATCH --time=00:30:00
#SBATCH --mail-user=am503@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --threads-per-core=1
#SBATCH --mem=250GB
#SBATCH --export=ALL

source /scratch/am503/miniforge3/etc/profile.d/conda.sh
conda activate seqscreen


nextflow run /scratch/am503/SOMAteM/subworkflows/seqscreen_v2.nf --input_dir /scratch/am503/SOMAteM/examples/data/fasta \
    --databases /scratch/am503/SeqScreenDB_23.4 \
    --working /scratch/am503/SOMAteM/examples/data/v3_seqscreen_workdir \
    --mode fast \
    --workflows_dir /scratch/am503/SOMAteM/subworkflows/local/ss2_subworkflows \
    --bin_dir /scratch/am503/SOMAteM/bin/seqscreen \
    --module_dir /scratch/am503/SOMAteM/bin/seqscreen/ss2_modules \
    --mode fast \
    --threads 96 \
    -c /scratch/am503/SOMAteM/conf/seqscreen.config \
    -with-conda \
    --sequencing_type long_read \
    -resume