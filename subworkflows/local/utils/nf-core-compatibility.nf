// Vendored-module compatibility processes to minimize boilerplate while testing

// Convert a filepath channel into a tuple including meta (id = filename, single_end = true)
// Supports: 
// - single files : makes [map, file]
// - directories with glob patterns : makes multiple [map, file] within a channel (?)
// - wildcard patterns (which get combined into a single tuple) : makes [map, multiple files]
workflow convert_to_nfcore_tuple {
    
    take:
    reads // string: path to reads

    main:

    // read single file
    is_multi_file = reads.endsWith("/")
    // is_single_fastaq = reads.endsWith(".fastq.gz") || reads.endsWith(".fastq") || reads.endsWith(".fa") || reads.endsWith(".fasta")
    combine_multiple_files = reads.contains("*")
    if (is_multi_file) {
        // read multiple files from directory and return a channel with multiple streams?
        tuple_out = channel.fromPath("${reads}/*.fastq.gz")
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by vendored modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    }
    if (combine_multiple_files) {
        // combine multiple files into a single tuple
        tuple_out = channel.fromPath(reads).collect()
            .map { files ->
                    def meta = [:] // Use dummy values; meta is required by vendored modules
                    meta.id = files.collect{f -> f.simpleName}.join("_")
                    meta.single_end = true
                    return [meta, files] }
        }
    else {
        tuple_out = channel.fromPath(reads)
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by vendored modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    }
    emit:
    tuple_out // tuple: [ meta, reads] of channels
}
