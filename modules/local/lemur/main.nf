// nf-core module file template copied from gms-16S/emu module
// refer there for more comments and stuff I deleted to maintain brevity..


process lemur {
    conda "bioconda::lemur"
    
    input:
    path reads
    path database_dir
    path taxonomy
    val rank

    output:
    path output_dir

    script:
    output_dir = "${reads.baseName}-lemur-output"
    """
    lemur -i $reads \
      -o $output_dir \
      -d $database_dir \
      --tax-path $taxonomy \
      -r $rank
    """
}

