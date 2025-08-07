#!/bin/bash
# Download Sylph database name from arguments
# Usage: ./download-db.sh <database_name>
# Reading: https://sylph-docs.github.io/pre%E2%80%90built-databases/

# If no arguments are provided, print usage and exit
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <database_name> such as gtdb-r226-c200-dbv1.syldb
    or smaller v0.3-c1000-gtdb-r214.syldb
    Read more at: https://sylph-docs.github.io/pre%E2%80%90built-databases/"
    exit 1
fi

# check if the file exists
if [ -f "$1" ]; then
    echo "File $1 already exists."
    exit 1
fi


cd ../../../databases # go to the databases directory

db_url="http://faust.compbio.cs.cmu.edu/sylph-stuff"

wget "${db_url}/$1"

