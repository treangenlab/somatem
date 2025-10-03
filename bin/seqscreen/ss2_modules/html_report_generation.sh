#!/bin/bash
set -e

usage() {
    echo;
    echo "Usage: $0 --report=/Path/to/seqscreen_report.txt --fasta=/Path/to/infile.fasta --blastx=/Path/to/blastx.xml --out=/Path/to/output.zip"
    echo "  --report      Path to seqscreen tsv report file"
    echo "  --fasta       Path to FASTA file used by seqscreen"
    echo "  --blastx      Path to BLASTx XML output file"
    echo "  --out         Output zip file"
    echo "  --mode        Sensitive/Default Mode"
    echo "  --rflag       Information about the enabled flags"
    echo "  --version     Version name (vX.X.X)"
    echo "  --funsocs     OPTIONAL Path to funsocs description txt file"
    echo "                defaults to scripts/html_report_generation/data/funsocs.txt"
    echo "  --gonames    OPTIONAL Path to go_names tsv file"
    echo "                defaults to /scratch1/external_databases/go_names.txt"
    echo "  --gonetwork  OPTIONAL Path to go_network tsv file"
    echo "  -h, --help  Print this help message"; echo;
    exit 1;
}

## ========================
## Globals

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINES="${DIR}"
SCRIPTS="${PIPELINES}/../scripts"
THREADS=1

## ========================
## Process input parameters
if [ $# -gt 8 ] || [ $# -lt 4 ]
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
    --mode=?*)
        MODE=${1#*=};;
    --mode|mode=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --version=?*)
        VERSION=${1#*=};;
    --version|version=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --report=?*)
        REPORT=${1#*=};;
    --report|report=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --rflag=?*)
	RFLAG=${1#*=};;
    --blastx=?*)
        BLASTX=${1#*=};;
    --blastx|blastx=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --funsocs=?*)
        FUNSOCS=${1#*=};;
    --gonames=?*)
        GONAMES=${1#*=};;
    --gonetwork=?*)
        GONETWORK=${1#*=};;
    --gonetwork|gonetwork=)
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

## Set defaults for optional arguments
if [ -z "$GONAMES" ]
then
    GONAMES="${SCRIPTS}/html_report_generation/data/go_names.txt"
fi

## Check if all required arguments are set
## Need to do that because there's optional arguments now
if [ -z "$FASTA" ] || [ -z "$REPORT" ] || [ -z "$BLASTX" ] || [ -z "$OUT" ]
then
    usage
fi
if [ -z "$FUNSOCS" ]
then
    FUNSOCS="${SCRIPTS}/html_report_generation/data/funsocs_description.txt"
fi

##===
## Logging
LOG="$( cd "$( dirname "${OUT}" )" && pwd )"
LOGTIME="${LOG}/html_report_generation.log"

## Run the the python code for HTML report generation
command time -v -o "${LOGTIME}" -a python3 \
    "${SCRIPTS}/html_report_generation/generateHtmlReport.py" \
    -f "${FASTA}" \
    -r "${REPORT}" \
    -b "${BLASTX}" \
    -o "${OUT}" \
    --rflag "${RFLAG}" \
    --version "${VERSION}" \
    -d "${SCRIPTS}/html_report_generation/libs/" \
    -t "${SCRIPTS}/html_report_generation/data/template.html" \
    -g "${GONAMES}" \
    -n "${GONETWORK}" \
    --funsocs "${FUNSOCS}" \
    --mode "${MODE}" \
    --go_template "${SCRIPTS}/html_report_generation/data/go_template.html" > "${LOGTIME}" 2>&1
status=$?

exit $status
