// nf-core compatibility processes to minimize boilerplate while testing

// Convert a filepath channel into a tuple including meta (id = filename, single_end = true)
workflow convert_to_nfcore_tuple {
    
    take:
    reads // string: path to reads

    main:

    // read single file
    is_single_file = reads.endsWith(".fastq.gz") || reads.endsWith(".fastq") || reads.endsWith(".fa") || reads.endsWith(".fasta")
    if (is_single_file) {
        tuple_out = Channel.fromPath(reads)
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    } else {

        // read multiple files from directory
        tuple_out = Channel.fromPath("${reads}/*.fastq.gz")
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
    }
    emit:
    tuple_out // tuple: [ meta, reads] of channels
}