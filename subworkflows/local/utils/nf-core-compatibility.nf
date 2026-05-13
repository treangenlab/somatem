// nf-core compatibility processes to minimize boilerplate while testing

// Convert a filepath channel into a tuple including meta (id = filename, single_end = true)
// Supports: single files, directories with glob patterns, and wildcard patterns (which get combined into a single tuple)
workflow convert_to_nfcore_tuple {
    
    take:
    reads // string: path to reads

    main:

    // read single file
    is_multi_file = reads.endsWith("/")
    is_single_file = reads.endsWith(".fastq.gz") || reads.endsWith(".fastq") || reads.endsWith(".fa") || reads.endsWith(".fasta")
    combine_multiple_files = reads.contains("*")
    if (is_multi_file) {
        // read multiple files from directory
        tuple_out = channel.fromPath("${reads}/*.fastq.gz")
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    }
    if (is_single_file) {
        tuple_out = channel.fromPath(reads)
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    }
    if (combine_multiple_files) {
        // combine multiple files into a single tuple
        tuple_out = channel.fromPath(reads).collect()
            .map { files ->
                    def meta = [:] // Use dummy values; meta is required by nf-core modules
                    meta.id = files.collect{f -> f.simpleName}.join("_")
                    meta.single_end = true
                    return [meta, files] }
        }
    else {
        // TODO: handle other cases
        error("Unsupported file format: ${reads}. Required: single file or directory with *.fastq.gz or *.fastq files")
    }
    emit:
    tuple_out // tuple: [ meta, reads] of channels
}