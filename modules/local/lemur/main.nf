// nf-core module file template copied from gms-16S/emu module
// refer there for more comments and stuff I deleted to maintain brevity..

reads = "examples/lemur/example-data/lemur_example.fastq"
database_dir = "examples/lemur/example-db"
taxonomy = "examples/lemur/example-db/taxonomy.tsv"
rank = "species"

output_dir = "examples/lemur/example-output"

process LEMUR {
    conda "bioconda::lemur"
    
    input:
    path(reads)

    output:
    path(output_dir)

    script:
    """
    lemur -i $reads \
      -o $output_dir \
      -d $database_dir \
      --tax-path $taxonomy \
      -r $rank
    """
}

