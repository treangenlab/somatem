// custom somatem module
process BAKTA_BAKTA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bakta=1.11.4"
    
    publishDir "${params.output_dir}/annotation/${meta.sample_id}/${meta.id}", mode: 'copy', pattern: "*.{embl,faa,ffn,fna,gbff,gff,tsv,txt,json,png,svg}"

    input:
    tuple val(meta), path(fasta)
    path(db)
    path(proteins)
    path(prodigal_tf)

    output:
    tuple val(meta), path("*.embl")           , emit: embl
    tuple val(meta), path("*.faa")            , emit: faa
    tuple val(meta), path("*.ffn")            , emit: ffn
    tuple val(meta), path("*.fna")            , emit: fna
    tuple val(meta), path("*.gbff")           , emit: gbff
    tuple val(meta), path("*.gff")            , emit: gff
    tuple val(meta), path("*.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("*.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("*.tsv")            , emit: tsv
    tuple val(meta), path("*.txt")            , emit: txt
    tuple val(meta), path("*.inference.tsv")  , emit: inference_tsv
    tuple val(meta), path("*.png")            , emit: png
    tuple val(meta), path("*.svg")            , emit: svg
    tuple val(meta), path("*.json")           , emit: json
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Use the original filename without extension as prefix
    def fasta_name = fasta.name
    def prefix = fasta_name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
    
    """
    # Check if FASTA file is valid and not empty
    if [ ! -s "${fasta}" ]; then
        echo "ERROR: Input FASTA file is empty or does not exist"
        exit 1
    fi
    
    # Check if FASTA has sequences
    seq_count=\$(grep -c '^>' ${fasta} || echo 0)
    if [ "\$seq_count" -eq 0 ]; then
        echo "ERROR: No sequences found in FASTA file"
        exit 1
    fi
    
    # Run Bakta annotation
    bakta \\
        --db ${db} \\
        --verbose \\
        --meta \\
        --output . \\
        --force \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fasta}
    
    bakta_exit_code=\$?
    
    if [ \$bakta_exit_code -ne 0 ]; then
        echo "ERROR: Bakta failed with exit code \$bakta_exit_code"
        exit \$bakta_exit_code
    fi
    
    # Define all expected output files
    expected_files=(
        "${prefix}.embl"
        "${prefix}.faa" 
        "${prefix}.ffn"
        "${prefix}.fna"
        "${prefix}.gbff"
        "${prefix}.gff"
        "${prefix}.hypotheticals.tsv"
        "${prefix}.hypotheticals.faa"
        "${prefix}.tsv"
        "${prefix}.txt"
        "${prefix}.inference.tsv"
        "${prefix}.json"
    )
    
    # Ensure all expected output files exist
    for file in "\${expected_files[@]}"; do
        if [ ! -f "\$file" ]; then
            echo "WARNING: Expected file \$file not found, creating empty file"
            touch "\$file"
        fi
    done

    # Generate visualization plots using bakta_plot
    if [ -f "${prefix}.json" ] && [ -s "${prefix}.json" ]; then
        bakta_plot \\
            ${prefix}.json \\
            --dpi 300 \\
            -o ${prefix}_plots 2>/dev/null || {
            echo "WARNING: bakta_plot failed, creating empty plot files"
            touch "${prefix}.png"
            touch "${prefix}.svg"
        }
        
        # Move the generated plots to the expected output names
        if [ -f "${prefix}_plots/${prefix}.png" ]; then
            mv "${prefix}_plots/${prefix}.png" "${prefix}.png"
        else
            touch "${prefix}.png"
        fi
        
        if [ -f "${prefix}_plots/${prefix}.svg" ]; then
            mv "${prefix}_plots/${prefix}.svg" "${prefix}.svg"
        else
            touch "${prefix}.svg"
        fi
        
        # Clean up the plots directory
        rm -rf "${prefix}_plots" 2>/dev/null || true
    else
        touch "${prefix}.png"
        touch "${prefix}.svg"
    fi
    
    # Final check - ensure all output files exist
    all_outputs="${prefix}.embl ${prefix}.faa ${prefix}.ffn ${prefix}.fna ${prefix}.gbff ${prefix}.gff ${prefix}.hypotheticals.tsv ${prefix}.hypotheticals.faa ${prefix}.tsv ${prefix}.txt ${prefix}.inference.tsv ${prefix}.png ${prefix}.svg ${prefix}.json"
    
    for file in \$all_outputs; do
        if [ ! -f "\$file" ]; then
            touch "\$file"  # Create it anyway to prevent pipeline failure
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version 2>&1) | sed 's/^.*bakta //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def fasta_name = fasta.name
    def prefix = fasta_name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
    """
    touch ${prefix}.embl
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.fna
    touch ${prefix}.gbff
    touch ${prefix}.gff
    touch ${prefix}.hypotheticals.tsv
    touch ${prefix}.hypotheticals.faa
    touch ${prefix}.tsv
    touch ${prefix}.txt
    touch ${prefix}.inference.tsv
    touch ${prefix}.png
    touch ${prefix}.svg
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version 2>&1) | sed 's/^.*bakta //; s/ .*\$//')
    END_VERSIONS
    """
}