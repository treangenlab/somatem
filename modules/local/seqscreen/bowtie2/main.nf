process SEQSCREEN_BOWTIE2 {
    tag "$meta.id:$target"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/bowtie2/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    val target

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_path = target == 'vfdb' ? "${db}/bowtie2/vfdb/vfdb" : "${db}/bowtie2/bsat_ccl/blacklist.seqs.nt"
    def out = target == 'vfdb' ? "blacklist_vfdb.sam" : "blacklist_bsat.sam"
    """
    bowtie2 \\
        --threads ${task.cpus} \\
        --sensitive \\
        -f \\
        --no-head \\
        --no-unal \\
        -x ${db_path} \\
        -U ${fasta} \\
        -S ${out} \\
        > bowtie2.log 2>&1 || {
            status=\$?
            cat bowtie2.log >&2
            # Recover the child signal when the wrapper masks a retryable failure.
            if [[ \$status == 1 ]]; then
                if grep -Fq 'bowtie2-align died with signal 9 (KILL)' bowtie2.log; then
                    exit 137
                fi
                if grep -Fq 'std::bad_alloc' bowtie2.log && grep -Fq 'bowtie2-align died with signal 6 (ABRT)' bowtie2.log; then
                    exit 134
                fi
            fi
            exit "\$status"
        }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: "\$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*version //')"
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_${target}.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: "stub"
    END_VERSIONS
    """
}
