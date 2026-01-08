/*
 * Parse MAGnet cluster representative CSV to extract reference IDs
 * Outputs a list of reference IDs that should be used for assembly
 */

process PARSE_MAGNET_REFS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(magnet_csv)

    output:
    tuple val(meta), path("selected_refs.txt"), emit: ref_ids
    path "versions.yml"                       , emit: versions

    script:
    """
    #!/usr/bin/env python3
    import csv
    import sys
    
    # Parse MAGnet cluster_representative.csv
    # Extract Assembly Accession IDs that are marked as "Present"
    # MAGnet output columns include:
    # - 'Assembly Accession ID': The genome accession (e.g., GCF_028867355.1)
    # - 'Presence/Absence': Status ('Present', 'Absent', 'Genus Present')
    # - 'Cluster Representative': Boolean if this genome is the cluster rep
    
    selected_refs = []
    
    try:
        with open('${magnet_csv}', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Extract references marked as Present
                presence = row.get('Presence/Absence', '')
                accession = row.get('Assembly Accession ID', '')
                
                # Include Present genomes (and optionally Genus Present)
                if presence == 'Present' and accession:
                    selected_refs.append(accession)
                    print(f"Selected: {accession} (Status: {presence})", file=sys.stderr)
                elif presence == 'Genus Present' and accession:
                    # Optionally include Genus Present matches
                    selected_refs.append(accession)
                    print(f"Selected: {accession} (Status: {presence})", file=sys.stderr)
                elif presence == 'Absent':
                    print(f"Skipped: {accession} (Status: {presence})", file=sys.stderr)
                    
    except Exception as e:
        print(f"Error: Could not parse MAGnet CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write selected reference IDs to output
    print(f"Total selected references: {len(selected_refs)}", file=sys.stderr)
    
    with open('selected_refs.txt', 'w') as out:
        for ref_id in selected_refs:
            out.write(f"{ref_id}\\n")
    
    # If no references were selected, this is a warning but not an error
    # The downstream process will need to handle empty files appropriately
    if len(selected_refs) == 0:
        print("WARNING: No references marked as Present", file=sys.stderr)
    
    # Version output
    with open('versions.yml', 'w') as v:
        v.write('"${task.process}":\\n')
        v.write('    python: "' + sys.version.split()[0] + '"\\n')
    """

    stub:
    """
    echo "ref1" > selected_refs.txt
    echo "ref2" >> selected_refs.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
