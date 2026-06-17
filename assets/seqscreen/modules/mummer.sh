#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --fasta=/Path/to/infile.fasta --database=/Path/to/NUC_DB --out=/Path/to/output.tab"
    echo "  --fasta     Path to the input NT FASTA file"
    echo "  --out       Path to the output BTAB file"
    echo "  --database  Path to the BLASTn database index base"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES=${DIR}

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
    --out=?*)
	    OUT=${1#*=};;
    --out|out=)
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
LOGTIME="${LOG}/mummer.log"

## Run MUMmer
command time -v -o ${LOGTIME} -a mummer -maxmatch -l 4 \
    "${FASTA}" \
    "${DATABASE}" \
    > "${OUT}" 2> "${LOGTIME}"
status=$?

exit $status
