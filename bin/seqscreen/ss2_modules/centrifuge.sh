#!/bin/bash
set -e

usage() {

    echo; echo "Usage: $0 -x=/Path/to/CentrifugeIndex -f=<toggle fasta> -U=<Comma separated reads to classify> 
    											-S=/Path/to/output.tsv  --threads=num_threads"
    echo "  -x             Path to the Centrifuge Index created"
    echo "  -f/-q          Read files are FASTA/FASTQ respectively"
    echo "  -U             Comma separated read files (reads1.fq,reads2.fq)"
    echo "  -S             Path to the output file (Default=stdout)"
    echo "  --threads      Number of threads to use (Default=1)"
    echo "  -h, --help     Print this help message out"; echo;
    exit 1;
}


## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES="${DIR}"
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


##============
## Logging
LOG="$( cd "$( dirname "${OUT}" )" && pwd )"
LOGTIME="${LOG}/centrifuge.log"

## Check number of CPUs on this machine and adjust if needed

## Run Centrifuge
command time -v -o "${LOGTIME}" -a centrifuge -x "${DATABASE}" \
	-f \
	-U "${FASTA}" \
	-S "${OUT}" \
	--threads "${THREADS}" \
        --seed 1 \
     > "${LOGTIME}" 2>&1
status=$?
  
exit $status

