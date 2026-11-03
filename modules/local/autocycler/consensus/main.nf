process AUTOCYCLER_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler/autocycler_out" },
        mode: params.publish_dir_mode,
        pattern: "autocycler_out/**"
    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/assembly" },
        mode: params.publish_dir_mode,
        pattern: "*.autocycler.{fasta,gfa}"
    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler" },
        mode: params.publish_dir_mode,
        pattern: "*.autocycler.{log,warning.txt}"

    input:
    tuple val(meta), path(assemblies), path(genome_size_file)

    output:
    tuple val(meta), path("*.autocycler.fasta"), emit: assembly, optional: true
    tuple val(meta), path("*.autocycler.gfa")  , emit: graph, optional: true
    tuple val(meta), path(assemblies), path("*.autocycler.warning.txt"), emit: fallback, optional: true
    tuple val(meta), path("*.autocycler.log"), emit: log
    tuple val(meta), path("autocycler_out")    , emit: outdir
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    mkdir -p assemblies
    for assembly in ${assemblies}; do
        cp "\${assembly}" assemblies/
    done

    # Run in a separate shell so errexit still applies inside this conditional.
    # Keep partial diagnostics even when compress/cluster/resolve/combine fails.
    if bash -euo pipefail > ${prefix}.autocycler.log 2>&1 <<'RUN_AUTOCYCLER'
    autocycler compress \\
        -i assemblies \\
        -a autocycler_out \\
        --max_contigs ${params.autocycler_max_contigs}

    autocycler cluster \\
        -a autocycler_out \\
        --max_contigs ${params.autocycler_max_contigs}

    shopt -s nullglob
    clusters=(autocycler_out/clustering/qc_pass/cluster_*)
    if [ \${#clusters[@]} -eq 0 ]; then
        echo "No Autocycler QC-pass clusters were produced for ${prefix}" >&2
        exit 1
    fi

    for cluster in "\${clusters[@]}"; do
        autocycler trim -c "\${cluster}"
        autocycler resolve -c "\${cluster}"
    done

    autocycler combine \\
        -a autocycler_out \\
        -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

    # A surviving small cluster does not establish that the genome was recovered.
    awk -v expected="\$(cat '${genome_size_file}')" \\
        -v minimum=${params.autocycler_min_genome_fraction} '
        BEGIN {
            expected = tolower(expected)
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", expected)
            scale = expected ~ /k\$/ ? 1000 : expected ~ /m\$/ ? 1000000 : expected ~ /g\$/ ? 1000000000 : 1
            sub(/[kmg]\$/, "", expected)
            if (expected !~ /^[+]?([0-9]+([.][0-9]*)?|[.][0-9]+)(e[+-]?[0-9]+)?\$/ || expected * scale < 0.5) {
                print "Invalid genome size for consensus validation: " expected > "/dev/stderr"
                invalid = 1
                exit 1
            }
            expected = int(expected * scale + 0.5)
        }
        !/^>/ { gsub(/[[:space:]]/, ""); bases += length(\$0) }
        END {
            if (invalid) exit 1
            printf "Consensus size: %.0f bp; expected: %.0f bp; recovered: %.2f%%\\n", bases, expected, 100 * bases / expected > "/dev/stderr"
            if (bases == 0 || bases < expected * minimum) {
                printf "Incomplete Autocycler consensus: below %.2f%% of the expected genome size.\\n", 100 * minimum > "/dev/stderr"
                print "Inspect input_assemblies.yaml and clustering/{qc_pass,qc_fail} in autocycler_out in this task work directory; a chromosome cluster may have failed QC." > "/dev/stderr"
                exit 1
            }
        }
    ' autocycler_out/consensus_assembly.fasta

    test -s autocycler_out/consensus_assembly.gfa
    RUN_AUTOCYCLER
    then
        cp autocycler_out/consensus_assembly.fasta ${prefix}.autocycler.fasta
        cp autocycler_out/consensus_assembly.gfa ${prefix}.autocycler.gfa
    else
        echo "Sample ${prefix}: Autocycler failed or produced an incomplete consensus. Selecting an input assembly using CheckM2 completeness - 5 * contamination. See ${prefix}.autocycler.log and autocycler_out for diagnostics." > ${prefix}.autocycler.warning.txt
        cat ${prefix}.autocycler.warning.txt >&2
    fi
    mkdir -p autocycler_out
    cat ${prefix}.autocycler.log >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | sed 's/^autocycler //')
    END_VERSIONS
    """
}
