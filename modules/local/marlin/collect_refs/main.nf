/*
 * Collect reference genomes downloaded by MAGnet
 * MAGnet downloads references to reference_genomes/ directory
 * This process collects those FASTA files for the selected references
 */

process COLLECT_MAGNET_REFS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(ref_ids), path(magnet_dir)

    output:
    tuple val(meta), path("selected_refs/*.fasta"), emit: refs
    path "versions.yml"                            , emit: versions

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    mkdir -p selected_refs
    
    # Read the list of selected reference IDs
    if [ ! -s ${ref_ids} ]; then
        echo "ERROR: No references selected by MAGnet" >&2
        exit 1
    fi
    
    # Copy selected reference FASTA files from MAGnet output
    # MAGnet stores references in reference_genomes/ subdirectory
    ref_dir=""
    if [ -d "${magnet_dir}/reference_genomes" ]; then
        ref_dir="${magnet_dir}/reference_genomes"
    elif [ -d "reference_genomes" ]; then
        ref_dir="reference_genomes"
    else
        echo "ERROR: Cannot find reference_genomes directory" >&2
        exit 1
    fi
    
    count=0
    while IFS= read -r accession; do
        if [ -n "\$accession" ]; then
            ref_file="\${ref_dir}/\${accession}.fasta"
            if [ -f "\$ref_file" ]; then
                cp "\$ref_file" "selected_refs/\${accession}.fasta"
                count=\$((count + 1))
                echo "Collected: \${accession}.fasta" >&2
            else
                echo "WARNING: Reference file not found: \$ref_file" >&2
            fi
        fi
    done < ${ref_ids}
    
    echo "Total references collected: \$count" >&2
    
    if [ \$count -eq 0 ]; then
        echo "ERROR: No reference files could be collected" >&2
        exit 1
    fi
    
    # Version output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p selected_refs
    echo ">stub_ref" > selected_refs/stub.fasta
    echo "ACGTACGT" >> selected_refs/stub.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: stub
    END_VERSIONS
    """
}
