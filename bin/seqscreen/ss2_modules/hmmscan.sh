#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --fasta=/Path/to/infile.fasta --database=/Path/to/FASTA_DB --out=/Path/to/output.txt --evalue=1e-5 [--threads=4]"
    echo "  --fasta     Path to the input peptide FASTA file"
    echo "  --out       Path to the output tabular file"
    echo "  --database  Path to the FASTA database file"
    echo "  --evalue    The E value cut-off for the HMMER analysis (Default = 1e-3)"
    echo "  --threads   Number of threads to use (Default=1)"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES="${DIR}"
SCRIPTS="${PIPELINES}/../scripts"
THREADS=1
EVALUE="1e-3"

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
LOGTIME="${LOG}/phmmer.log"

## This isn't the pHMMER binary, this is a wrapper script that runs
## an embarassingly parallel instance of hmmscan
command time -v -o "${LOGTIME}" -a "${SCRIPTS}/para_hmmscan.pl" \
    -q "${FASTA}" \
    -d "${DATABASE}" \
    -o "${OUT}" \
    -e "${EVALUE}" \
    -t "${THREADS}" > "${LOGTIME}" 2>&1
status=$?

exit $status
