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
        pattern: "*.autocycler.{log,lengths.tsv}"

    input:
    tuple val(meta), path(assemblies), path(genome_size_file)

    output:
    tuple val(meta), path("*.autocycler.fasta"), emit: assembly
    tuple val(meta), path("*.autocycler.gfa")  , emit: graph
    tuple val(meta), path("*.autocycler.log"), emit: log
    tuple val(meta), path("*.autocycler.lengths.tsv"), emit: lengths
    tuple val(meta), path("autocycler_out")    , emit: outdir
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def assembly_count = assemblies instanceof Collection ? assemblies.size() : 1
    def min_assemblies_arg = params.autocycler_min_assemblies == null ? '' : "--min_assemblies ${params.autocycler_min_assemblies}"
    """
    # Stop on any failed step while streaming and retaining diagnostics.
    set -o pipefail
    bash -euo pipefail <<'RUN_AUTOCYCLER' 2>&1 | tee ${prefix}.autocycler.log
    expected_genome_size=\$(cat '${genome_size_file}')
    # Apply the same whole-assembly length rule before clustering and after consensus.
    check_length() {
        awk -v expected="\$expected_genome_size" \\
            -v minimum=${params.autocycler_min_genome_fraction} \\
            -v stage="\$1" -v assembly="\$2" '
            BEGIN {
                expected = tolower(expected)
                gsub(/^[[:space:]]+|[[:space:]]+\$/, "", expected)
                scale = expected ~ /k\$/ ? 1000 : expected ~ /m\$/ ? 1000000 : expected ~ /g\$/ ? 1000000000 : 1
                sub(/[kmg]\$/, "", expected)
                if (expected !~ /^[+]?([0-9]+([.][0-9]*)?|[.][0-9]+)(e[+-]?[0-9]+)?\$/ || expected * scale < 0.5) {
                    print "Invalid genome size for length validation: " expected > "/dev/stderr"
                    invalid = 1
                    exit 2
                }
                expected = int(expected * scale + 0.5)
            }
            !/^>/ { gsub(/[[:space:]]/, ""); bases += length(\$0) }
            END {
                if (invalid) exit 2
                passes = bases > 0 && bases >= expected * minimum
                printf "%s\\t%s\\t%.0f\\t%.0f\\t%.6f\\t%s\\n", stage, assembly, bases, expected, bases / expected, passes ? "True" : "False"
                printf "%s %s: %.0f bp; expected: %.0f bp; recovered: %.2f%%; passes length filter: %s\\n", stage, assembly, bases, expected, 100 * bases / expected, passes ? "yes" : "no" > "/dev/stderr"
                exit (passes ? 0 : 1)
            }
        ' "\$2"
    }

    mkdir -p assemblies
    printf 'Stage\\tAssembly\\tLength_bp\\tExpected_bp\\tGenome_fraction\\tPasses_length_filter\\n' > ${prefix}.autocycler.lengths.tsv
    for assembly in ${assemblies}; do
        if check_length candidate "\$assembly" >> ${prefix}.autocycler.lengths.tsv; then
            cp "\$assembly" assemblies/
        else
            check_status=\$?
            if [ "\$check_status" -ne 1 ]; then exit "\$check_status"; fi
        fi
    done
    shopt -s nullglob
    selected=(assemblies/*.fasta)
    echo "Autocycler inputs: \${#selected[@]} of ${assembly_count} candidates pass the length filter."
    if [ \${#selected[@]} -eq 0 ]; then
        echo "No candidate assemblies meet the minimum genome fraction for ${prefix}. Inspect ${prefix}.autocycler.lengths.tsv and the supplied or estimated genome size." >&2
        exit 1
    fi

    autocycler compress \\
        -i assemblies \\
        -a autocycler_out \\
        --max_contigs ${params.autocycler_max_contigs} \\
        --threads ${task.cpus}

    autocycler cluster \\
        -a autocycler_out \\
        --max_contigs ${params.autocycler_max_contigs} \\
        ${min_assemblies_arg}

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

    if ! check_length consensus autocycler_out/consensus_assembly.fasta >> ${prefix}.autocycler.lengths.tsv; then
        echo "Incomplete Autocycler consensus: below the minimum genome fraction. Inspect input_assemblies.yaml and clustering/{qc_pass,qc_fail} in autocycler_out in this task work directory." >&2
        exit 1
    fi

    test -s autocycler_out/consensus_assembly.gfa
    RUN_AUTOCYCLER

    cp autocycler_out/consensus_assembly.fasta ${prefix}.autocycler.fasta
    cp autocycler_out/consensus_assembly.gfa ${prefix}.autocycler.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | sed 's/^autocycler //')
    END_VERSIONS
    """
}
