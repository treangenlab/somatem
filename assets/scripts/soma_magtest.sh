#! /bin/bash

#SBATCH --job-name=somatem
#SBATCH --partition=debug
#SBATCH --account=commons
#SBATCH --time=00:30:00
#SBATCH --mail-user=am503@rice.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --threads-per-core=1
#SBATCH --mem=128GB
#SBATCH --export=ALL

#source /scratch/am503/miniforge3/etc/profile.d/conda.sh
#conda activate somatem

nextflow run /scratch/am503/SOMAteM/subworkflows/local/assembly_mags_sub.nf \
        --input /scratch/am503/SOMAteM/examples/data/input4mags \
        --outdir /scratch/am503/SOMAteM/examples/data/mag_test8 \
        --threads 96 \
        --flye_mode nano-hq \
        --semibin_environment mouse_gut \
        --completeness_threshold 80 \
        --checkm2_db /scratch/am503/aviary_dbs/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
        --bakta_db /scratch/am503/aviary_dbs/bakta/db \
        --singlem_db /scratch/am503/aviary_dbs/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb \
        -c /scratch/am503/SOMAteM/conf/nots_mag.config \
        -with-report /scratch/am503/SOMAteM/examples/data/mag_test8/reports/mag_test_report.html \
        -with-trace /scratch/am503/SOMAteM/examples/data/mag_test8/reports/mag_test_trace.txt \
        -resume
