process AUTOCYCLER_FALLBACK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../../nf-core/checkm2/predict/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/assembly" },
        mode: params.publish_dir_mode,
        pattern: "*.fallback.fasta"
    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/quality/assembly_selection" },
        mode: params.publish_dir_mode,
        pattern: "*.tsv"

    input:
    tuple val(meta), path(assemblies), path(warning)
    path(db)

    output:
    tuple val(meta), path("*.fallback.fasta"), emit: assembly
    tuple val(meta), path("*.assembly_selection.tsv"), emit: selection
    tuple val(meta), path("*.candidate_quality.tsv"), emit: quality
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    checkm2 predict \\
        --input ${assemblies} \\
        --output-directory candidate_quality \\
        --threads ${task.cpus} \\
        --database_path ${db}

    cp candidate_quality/quality_report.tsv ${prefix}.candidate_quality.tsv
    python - '${prefix}' ${assemblies} <<'SELECT_ASSEMBLY'
    import csv
    import math
    from pathlib import Path
    import shutil
    import sys

    sample = sys.argv[1]
    assemblies = {Path(name).stem: Path(name) for name in sys.argv[2:]}
    if len(assemblies) != len(sys.argv[2:]):
        raise ValueError('Duplicate candidate assembly names')
    with open('candidate_quality/quality_report.tsv') as handle:
        rows = list(csv.DictReader(handle, delimiter='\\t'))
    if len(rows) != len(assemblies) or {row['Name'] for row in rows} != set(assemblies):
        raise ValueError('CheckM2 must report quality for every candidate assembly')
    for row in rows:
        completeness = float(row['Completeness'])
        contamination = float(row['Contamination'])
        if not (math.isfinite(completeness) and 0 <= completeness <= 100
                and math.isfinite(contamination) and contamination >= 0):
            raise ValueError(f'Invalid CheckM2 quality for {row["Name"]}')
        row['Completeness'] = completeness
        row['Contamination'] = contamination
        row['Score'] = completeness - 5 * contamination
    # Stable ties: higher completeness, lower contamination, then filename.
    rows.sort(key=lambda row: (-row['Score'], -row['Completeness'], row['Contamination'], row['Name']))
    best = rows[0]
    selected = assemblies[best['Name']]
    if not selected.is_file() or selected.stat().st_size == 0:
        raise ValueError(f'Selected assembly is missing or empty: {selected}')
    shutil.copyfile(selected, f'{sample}.fallback.fasta')
    with open(f'{sample}.assembly_selection.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\\t')
        writer.writerow(['Sample', 'Assembly', 'Completeness', 'Contamination', 'Score', 'Selected'])
        for row in rows:
            writer.writerow([sample, assemblies[row['Name']].name, row['Completeness'],
                             row['Contamination'], row['Score'], row is best])
    print(f'WARNING: Sample {sample}: using {selected.name} after Autocycler failure '
          f'(completeness={best["Completeness"]}%, contamination={best["Contamination"]}%, '
          f'score={best["Score"]}).', file=sys.stderr)
    SELECT_ASSEMBLY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
