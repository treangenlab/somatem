// nf-core compatibility processes to minimize boilerplate while testing

// Convert a filepath channel into a tuple including meta
workflow convert_to_nfcore_tuple {
    
    take:
    reads // string: path to reads

    main:

    tuple_out = Channel.fromPath(reads)
            .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = "test"
                meta.single_end = false
                return [meta, r] }

    emit:
    tuple_out // tuple: [ meta, reads] of channels
}