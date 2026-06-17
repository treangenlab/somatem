#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --fasta=/Paht/to/infile.fasta --btab=/Path/to/infile.btab --out=/Path/to/output_directory"
    echo "  --fasta     Path to the input FASTA file used in the BLAST"
    echo "  --btab      Path to the input tabular BLAST output file"
    echo "  --out       Path to the output file"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES="${DIR}"
SCRIPTS="${PIPELINES}/../scripts"
THREADS=1

## ========================
## Process the input parameters
if [ $# -gt 3 ] || [ $# -lt 3 ]
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
    --btab=?*)
	    BTAB=${1#*=};;
    --btab|btab=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --out=?*)
	    OUT=${1#*=};;
    --out|out=)
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
LOGTIME="${LOG}/outlier_detection.log"

command time -v -o "${LOGTIME}" -a "sort" \
     -r  \
     -k "1,1" \
     -k "12,12" \
    "${BTAB}" > "${BTAB}.srt" 
status=$?

## Outlier detection commands below...
command time -v -o "${LOGTIME}" -a "${SCRIPTS}/outlier_detection/score_blast.py" \
    -q "${FASTA}" \
    -b "${BTAB}.srt" \
    -o "${OUT}" > "${LOGTIME}" 2>&1
status=$?

## Second script to create a cleans the BTAB file of outliers
command time -v -o "${LOGTIME}" -a "${SCRIPTS}/outlier_detection/remove_outliers.pl" \
    --outliers "${OUT}" \
    --btab "${BTAB}" \
    --out "${BTAB}_outlier_clean.btab"
status=$?

exit $status
