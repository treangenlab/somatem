/*
 * Filter reference sequences based on MAGnet selection
 * Keeps only the sequences for selected reference IDs
 */

process FILTER_REFERENCE_SEQS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(ref_ids), path(all_refs)

    output:
    tuple val(meta), path("filtered_refs.fasta"), emit: filtered_refs
    path "versions.yml"                         , emit: versions

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import sys
    
    # Load selected reference IDs
    selected_ids = set()
    try:
        with open('${ref_ids}', 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    selected_ids.add(line)
    except Exception as e:
        print(f"Warning: Could not read reference IDs: {e}", file=sys.stderr)
    
    # Filter reference sequences
    filtered_count = 0
    with open('filtered_refs.fasta', 'w') as out:
        for record in SeqIO.parse('${all_refs}', 'fasta'):
            # Check if sequence ID matches any selected reference
            seq_id = record.id
            # Try matching with or without common prefixes
            match_id = seq_id
            for prefix in ['ref|', 'gi|', 'gb|', 'emb|', 'dbj|']:
                if match_id.startswith(prefix):
                    match_id = match_id.split('|')[1] if '|' in match_id else match_id
                    break
            
            # Check if this sequence should be kept
            if seq_id in selected_ids or match_id in selected_ids:
                SeqIO.write(record, out, 'fasta')
                filtered_count += 1
            else:
                # Also check if any selected ID is a substring of seq_id
                for sel_id in selected_ids:
                    if sel_id in seq_id or seq_id in sel_id:
                        SeqIO.write(record, out, 'fasta')
                        filtered_count += 1
                        break
    
    print(f"Filtered {filtered_count} reference sequences from MAGnet selection", file=sys.stderr)
    
    # If no sequences matched, copy all references as fallback
    if filtered_count == 0:
        print("Warning: No sequences matched MAGnet selection, using all references", file=sys.stderr)
        import shutil
        shutil.copy('${all_refs}', 'filtered_refs.fasta')
    
    # Version output
    with open('versions.yml', 'w') as v:
        v.write('"${task.process}":\\n')
        v.write('    python: "' + sys.version.split()[0] + '"\\n')
        v.write('    biopython: "' + __import__('Bio').__version__ + '"\\n')
    """

    stub:
    """
    cp ${all_refs} filtered_refs.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        biopython: stub
    END_VERSIONS
    """
}
