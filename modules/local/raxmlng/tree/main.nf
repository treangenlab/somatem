process RAXMLNG_TREE {
    tag "raxml-ng"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:2.0.1--h3f2fef4_0':
        'quay.io/biocontainers/raxml-ng:2.0.1--h3f2fef4_0' }"

    input:
    path alignment

    output:
    path "core_snps.raxml.bestTree", emit: tree
    path "raxmlng/*"               , emit: raxmlng_dir
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode = params.snp_phylo_raxmlng_bootstrap ? "--all" : "--search"
    def args = task.ext.args ?: ''
    """
    mkdir -p raxmlng

    python3 - <<'PY'
    from pathlib import Path
    records = {}
    name = None
    for line in Path("${alignment}").read_text().splitlines():
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line.strip().upper())
    seqs = {k: "".join(v) for k, v in records.items()}
    variable = False
    if seqs:
        length = len(next(iter(seqs.values())))
        for idx in range(length):
            if len({seq[idx] for seq in seqs.values()}) > 1:
                variable = True
                break
    if not variable:
        names = list(seqs) or ["reference"]
        Path("core_snps.raxml.bestTree").write_text("(" + ",".join(names) + ");\\n")
        Path("raxmlng/core_snps.raxml.log").write_text("No variable core SNPs; wrote a star tree.\\n")
    PY

    if [[ ! -s core_snps.raxml.bestTree ]]; then
        raxml-ng \\
            ${mode} \\
            --msa ${alignment} \\
            --model ${params.snp_phylo_raxmlng_model} \\
            --prefix raxmlng/core_snps \\
            --threads ${task.cpus} \\
            --tree pars{${params.snp_phylo_raxmlng_parsimony_trees}} \\
            ${args}
        cp raxmlng/core_snps.raxml.bestTree core_snps.raxml.bestTree
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxml-ng: \$(raxml-ng --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p raxmlng
    echo "(reference,sample);" > core_snps.raxml.bestTree
    touch raxmlng/core_snps.raxml.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxml-ng: "stub"
    END_VERSIONS
    """
}
