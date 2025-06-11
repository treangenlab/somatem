// nf-core module file template copied from gms-16S/emu module
// 
//  Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
//  Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process LEMUR {
    debug true
    tag "$meta.id"
    label 'process_high'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda "bioconda::lemur"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bioconductor-lemur%3A1.4.0--r44he5774e6_0':
    //     'quay.io/biocontainers/emu:3.5.1--hdfd78af_0' }"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*abundance.tsv")                     , emit: report
    tuple val(meta), path("*read-assignment-distributions.tsv") , emit: assignment_report, optional:true
    tuple val(meta), path("*.sam")                              , emit: samfile, optional:true
    tuple val(meta), path("*.fastq_unclassified_mapped.fasta")  , emit: unclassified_mapped_fa , optional:true
    tuple val(meta), path("*.fastq_unmapped.fasta")             , emit: unclassified_unmapped_fa , optional:true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    lemur -i $reads \
      -o example-output \
      -d examples/example-db \
      --tax-path examples/example-db/taxnomy.tsv \
      -r species
    """
}

