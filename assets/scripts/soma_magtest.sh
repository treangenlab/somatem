#! /bin/bash

#SBATCH --job-name=somatem
#SBATCH --partition=commons
#SBATCH --account=commons
#SBATCH --time=23:30:00
#SBATCH --mail-user=am503@rice.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --threads-per-core=1
#SBATCH --mem=128GB
#SBATCH --export=ALL

source /scratch/am503/miniforge3/etc/profile.d/conda.sh
conda activate seqscreen

#nextflow run /scratch/am503/SOMAteM/subworkflows/mags_new.nf \
#        --input_dir /scratch/am503/SOMAteM/examples/data/input4mags \
#        --output_dir /scratch/am503/SOMAteM/examples/data/mag_test_final \
#        --threads 96 \
#        --flye_mode nano-hq \
#        --semibin_environment mouse_gut \
#        --completeness_threshold 50 \
#        --checkm2_db /scratch/am503/aviary_dbs/checkm2/uniref100.KO.1.dmnd \
#        --bakta_db /scratch/am503/aviary_dbs/bakta/db \
#        --singlem_metapackage /scratch/am503/aviary_dbs/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb \
#        -c /scratch/am503/SOMAteM/conf/nots_mag.config \
#        -with-report /scratch/am503/SOMAteM/examples/data/mag_test3/reports/mag_test_report.html \
#        -with-trace /scratch/am503/SOMAteM/examples/data/mag_test3/reports/mag_test_trace.txt \
#        -resume

# refassembly.nf subworkflow

#nextflow run /scratch/am503/SOMAteM/subworkflows/refassembly.nf \
#        --input_dir /scratch/am503/SOMAteM/examples/data/input4mags \
#        --output_dir /scratch/am503/SOMAteM/examples/data/mag_test_refassem \
#        --threads 80 \
#        -c /scratch/am503/SOMAteM/conf/refassembly.config \
#        --hostile_index "human-t2t-hla" \
#        --lemur_db /scratch/am503/dbs/rv221bacarc-rv222fungi \
#        --taxonomy /scratch/am503/dbs/rv221bacarc-rv222fungi/taxonomy.tsv \
#        --platform ont \
#        --caller medaka \
#        --iterations 3

# Set input and output base directories
#SEQSCREEN_DB="/scratch/am503/SeqScreenDB_23.4"

nextflow run /scratch/am503/SOMAteM/subworkflows/seqscreen.nf --fasta /scratch/am503/SOMAteM/examples/data/input4mags/fasta/s4_final.fasta \
	--databases /scratch/am503/SeqScreenDB_23.4 \
        --working /scratch/am503/SOMAteM/examples/data/seqscreen_workdir \
        --mode fast \
        --workflows_dir /scratch/am503/SOMAteM/subworkflows/local/ss2_subworkflows \
        --bin_dir /scratch/am503/SOMAteM/bin/seqscreen \
        --module_dir /scratch/am503/SOMAteM/bin/seqscreen/ss2_modules \
        --mode fast \
        --threads 96 \
        -c /scratch/am503/SOMAteM/conf/seqscreen.config \
        -with-conda \
        -resume