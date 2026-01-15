#!/bin/bash

# Get database version and information for a given taxonomy
# Exit codes:
#   0 - Success
#   1 - Invalid arguments
#   2 - CSV file not found
#   3 - Taxonomy not found

set -euo pipefail

# Print usage and exit with error
usage() {
    echo "Usage: $0 <taxonomy> <csv_path>" >&2
    exit 1
}

# Check arguments
[ $# -eq 2 ] || usage

target_taxonomy="$1"
csv_path="$2"

# Check if CSV file exists
[ -f "$csv_path" ] || { echo "Error: CSV file not found: $csv_path" >&2; exit 2; }

# Process CSV file
result=$(awk -F, -v target="$target_taxonomy" '
BEGIN { OFS="\t" }
NR == 1 { next }  # Skip header
$1 == target { 
    gsub(/^[\"\s]+|[\"\s]+$/,"",$2); 
    comment = $3; 
    gsub(/^[\"\s]+|[\"\s]+$/,"",comment); 
    print $2, (comment ? comment : "No description available"); 
    exit 0 
}
END { if (!found) exit 3 }' "$csv_path" 2>/dev/null)

# Handle results
case $? in
    0)  # Found
        db_version="${result%%$'\t'*}"
        db_comment="${result#*$'\t'}"
        echo "$db_version"
        ;;
    3)  # Not found
        available=$(tail -n +2 "$csv_path" | cut -d, -f1 | tr -d '"' | xargs | tr ' ' ',')
        echo "Error: Unknown target_taxonomy: $target_taxonomy" >&2
        echo "Available taxonomies: $available" >&2
        exit 3
        ;;
    *)  # Other errors
        echo "Error: Failed to process CSV file" >&2
        exit 1
        ;;
esac
