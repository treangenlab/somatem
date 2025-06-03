#! /bin/bash

#SBATCH --time=04-00:00:00
#SBATCH --partition=defq
#SBATCH --mail-user=email@myemail.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256GB
#SBATCH --nodes=1
#SBATCH --job-name=somatem
#SBATCH --comment=somatem

module load mamba

nextflow run /path/to/SOMAteM/workflows/somatem_prep.nf --input_dir /path/to/SOMAteM/examples/data --output_dir /path/to/SOMAteM/examples/soma_prep_out --threads 12 --maxlength 30000 --minq 10 --minlen 250 --host_index 'human-t2t-hla' 

