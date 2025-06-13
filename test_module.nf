workflow test {
    main:
    import LEMUR from "modules/local/lemur/main.nf"

    LEMUR(reads="examples/lemur/example-data/lemur_example.fastq")
}