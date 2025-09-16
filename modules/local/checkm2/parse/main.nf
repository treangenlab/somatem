// custom somatem module
process CHECKM2_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(checkm2_tsv)
    tuple val(meta2), path(bins)

    output:
    tuple val(meta), path("bins_with_completeness.csv"), emit: completeness_map
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create Python script
    cat > parse_checkm2.py << 'EOF'
import pandas as pd
import csv
import os

# Read CheckM2 results
checkm2_df = pd.read_csv('${checkm2_tsv}', sep='\\t')

# Create a mapping of bin name to completeness
completeness_map = {}
for _, row in checkm2_df.iterrows():
    bin_name = row['Name']
    completeness = float(row['Completeness'])
    completeness_map[bin_name] = completeness
    print(f"Bin {bin_name}: {completeness}% complete")

# Write the completeness mapping to a CSV file
with open('bins_with_completeness.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['bin_name', 'completeness'])
    for bin_name, completeness in completeness_map.items():
        writer.writerow([bin_name, completeness])

print(f"Created completeness mapping for {len(completeness_map)} bins")
EOF

    # Run the Python script
    python3 parse_checkm2.py

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch bins_with_completeness.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
