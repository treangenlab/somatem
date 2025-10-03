#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --fasta=/Path/to/infile.fasta --database=/Path/to/BLASTp_DB --out=/Path/to/output.btab --evalue=1e-5 [--threads=4]"
    echo "  --fasta     Path to the input nucleotide FASTA file"
    echo "  --out       Path to the output tabular file"
    echo "  --database  Path to the BLASTx database index base"
    echo "  --evalue    The E value cut-off for the BLAST analysis"
    echo "  --threads   Number of threads to use (Default=1)"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES=${DIR}
SCRIPTS=${PIPELINES}/../scripts
THREADS=1

## ========================
## Process the input parameters
if [ $# -gt 5 ] || [ $# -lt 4 ]
then
    usage
fi

while true
do
    case $1 in
    --help|-h)
	    usage
	    exit;;
    --fasta=?*)
	    FASTA=${1#*=};;
    --fasta|fasta=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --out=?*)
	    OUT=${1#*=};;
    --out|out=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --threads=?*)
	    THREADS=${1#*=};;
    --threads|threads=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --database=?*)
            DATABASE=${1#*=};;
    --database|database=)
            echo "$0: missing argument for '$1' option"
            usage
            exit 1;;
    --evalue=?*)
	    EVALUE=${1#*=};;
    --evalue|evalue=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --)
	    shift
	    break;;
    -?*)
	    echo "$0: invalid option: $1"
	    usage
	    exit 1;;
    *)
	    break
    esac
    shift
done

##===
## Logging
LOG="$( cd "$( dirname "${OUT}" )" && pwd )"
LOGTIME="${LOG}/blastx-fast.log"

## Run the BLASTx and produce an ASN file
command time -v -o ${LOGTIME} -a blastx -query "${FASTA}" \
    -db "${DATABASE}" \
    -out "${OUT}.asn" \
    -evalue "${EVALUE}" \
    -num_threads "${THREADS}" \
    -max_target_seqs 500 \
    -task blastx \
    -seg no \
    -outfmt 11 > "${LOGTIME}" 2>&1
status=$?

#    -max_target_seqs 500 \
## Convert the ASN file to a tabular and XML output in parallel
command time -v -o ${LOGTIME} -a ${SCRIPTS}/blast_formatter_parallel.pl -a ${OUT}.asn -o ${OUT} -f 5,"\"6 std ppos qframe score salltitles\""

rm ${OUT}.asn ## Clean up the ASN file, it's large and we don't need it anymore.

exit $status
