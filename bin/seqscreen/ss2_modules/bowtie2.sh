#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --fasta=/Path/to/infile.fasta --database=/Path/to/BOWTIE_DB --out=/Path/to/output.sam [--threads=4]"
    echo "  --fasta     Path to the input NT FASTA file"
    echo "  --out       Path to the output SAM file"
    echo "  --database  Path to the Bowtie2 database index base"
    echo "  --threads   Number of threads the pipeline use (Default=1)"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES="${DIR}"
THREADS=1

## ========================
## Process the input parameters
# [ $HOUR -gt 7 -a $HOUR -lt 17 ]
if [ $# -gt 4 ] || [ $# -lt 3 ]
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

##===
## Logging
LOG="$( cd "$( dirname "${OUT}" )" && pwd )"
LOGTIME="${LOG}/bowtie2.log"

## Run Bowtie2
command time -v -o ${LOGTIME} -a bowtie2 --threads "${THREADS}" --sensitive -f --no-head --no-unal\
    -x "${DATABASE}" \
    -U "${FASTA}" \
    -S "${OUT}" > "${LOGTIME}" 2>&1
status=$?

exit $status
