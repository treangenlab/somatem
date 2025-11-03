process VCONTACT2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vcontact2=0.11.3"

    input:
    tuple val(meta), path(viral_contigs)

    output:
    tuple val(meta), path("vcontact2_output/genome_by_genome_overview.csv"), emit: overview
    tuple val(meta), path("vcontact2_output/c1.ntw"), emit: network
    tuple val(meta), path("vcontact2_output/viral_cluster_overview.csv"), emit: clusters
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # First run gene prediction
    prodigal -i ${viral_contigs} -a proteins.faa -d genes.fna -f gff -o genes.gff

    # Create gene to genome mapping
    python3 <<EOF
import pandas as pd
from Bio import SeqIO

# Create gene to genome mapping
gene_to_genome = []
for record in SeqIO.parse("${viral_contigs}", "fasta"):
    genome_id = record.id
    # Assume prodigal naming convention
    for i in range(1, 1000):  # arbitrary large number
        gene_id = f"{genome_id}_{i}"
        gene_to_genome.append([gene_id, genome_id])

df = pd.DataFrame(gene_to_genome, columns=['protein_id', 'contig_id'])
df.to_csv('gene_to_genome.csv', index=False)
EOF

    vcontact2 \\
        --raw-proteins proteins.faa \\
        --rel-mode 'Diamond' \\
        --proteins-fp gene_to_genome.csv \\
        --db 'ProkaryoticViralRefSeq211-Merged' \\
        --output-dir vcontact2_output \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact2: \$(vcontact2 --version 2>&1 | grep 'vConTACT2 v' | sed 's/vConTACT2 v//g')
    END_VERSIONS
    """
}