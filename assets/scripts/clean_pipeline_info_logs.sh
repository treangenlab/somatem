#!/bin/bash
# This script cleans up the results/pipeline_info/ directory by removing files
# that do not match any timestamps found in the `nextflow log` output.
# First you need to clean the results directory. 
# You can use `nextflow clean -but funny_name1 -but name2 .. -f` ; by retaining the desired runs.

# run this from the root of the nextflow project directory in the nextflow micromamba environment
# Usage: ./archive/clean_pipeline_info_logs.sh [-d|--dry-run] [-h|--help]

# error: some issue with date conversion and grep parsing initially
# date: invalid date ‘2025-10-02 14-38:52’

# Help message
show_help() {
    echo "Usage: $0 [-d|--dry-run] [-h|--help]"
    echo "Options:"
    echo "  -d, --dry-run    Show what would be deleted without actually removing files"
    echo "  -h, --help       Show this help message"
    exit 0
}

# Parse command line arguments
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Function to convert timestamp to epoch seconds
to_epoch() {
    local timestamp="$1"
    # Convert both formats: "2025-10-03 00:15:20" and "2025-10-03_00-15-20"
    timestamp=$(echo "$timestamp" | sed 's/_/ /g')
    # Replace '_' with ' ', and the last three '-' in the time portion with ':'
    local formatted=$(echo "$timestamp" | sed 's/_/ /' | sed -E 's/([0-9]{4}-[0-9]{2}-[0-9]{2}) ([0-9]{2})-([0-9]{2})-([0-9]{2})/\1 \2:\3:\4/')
    date -d "$formatted" +%s
}

# create test cases for this to_epoch function
# echo "Testing to_epoch function..."
# test_dates=("2025-10-03 00:15:20" "2025-10-03_00-15-20" "2025-12-31")
# for date in "${test_dates[@]}"; do

# Function to convert nextflow log timestamp to pipeline_info format
convert_timestamp() {
    # Convert from format like "2025-10-03 00:15:20" to "2025-10-03_00-15-20"
    echo "$1" | sed 's/ /_/g' | sed 's/:/-/g'
}

# Function to check if two timestamps are within 2 minutes of each other
timestamps_match() {
    local ts1="$1"
    local ts2="$2"
    local epoch1=$(to_epoch "$ts1")
    local epoch2=$(to_epoch "$ts2")
    local diff=$((epoch1 - epoch2))
    diff=${diff#-} # absolute value
    [ "$diff" -le 5 ] # 5 seconds
}

# Get list of all timestamps from `nextflow log` output
echo "Extracting timestamps from `nextflow log` output..."
nextflow log > archive/temp_nextflow_log
timestamps=$(grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\} [0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}' archive/temp_nextflow_log | sort -u)

if [ -z "$timestamps" ]; then
    echo "No timestamps found in .nextflow.log"
    exit 1
fi

# Store timestamps in array
readarray -t timestamp_array <<< "$timestamps"

# Debug
echo "Valid timestamps extracted:"
for ts in "${timestamp_array[@]}"; do
    echo "$(convert_timestamp "$ts")"
done

# Process files in pipeline_info directory
if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Checking files in results/pipeline_info/..."
else
    echo "Checking files in results/pipeline_info/..."
fi

for file in results/pipeline_info/*; do
    if [ -f "$file" ]; then
        # Extract timestamp from filename (assuming format like execution_trace_2025-10-03_00-15-20.txt)
        file_timestamp=$(echo "$file" | grep -o '[0-9]\{4\}[-_][0-9]\{2\}[-_][0-9]\{2\}_[0-9]\{2\}-[0-9]\{2\}-[0-9]\{2\}')
        #  2025-10-03_00-15-20 becomes :        '2025       -  10         -   03       _ 00       - 15       - 20'
        
        # debug
        echo "timestamp from file: $file_timestamp" # print the timestamp

        if [ ! -z "$file_timestamp" ]; then
            # Convert file timestamp format for comparison
            file_timestamp_std=$(echo "$file_timestamp" | sed 's/_/ /g' | sed 's/-\([0-9][0-9]\)$/:\1/g')
            
            # Check if this timestamp matches any in our valid list within 2 minutes
            matched=0
            for valid_ts in "${timestamp_array[@]}"; do
                if timestamps_match "$valid_ts" "$file_timestamp_std"; then
                    matched=1
                    break
                fi
            done

            if [ $matched -eq 1 ]; then
                if [ $DRY_RUN -eq 1 ]; then
                    echo "[DRY RUN] Would keep: $file (matched within 2 minutes)"
                else
                    echo "Keeping matched file: $file"
                fi
            else
                if [ $DRY_RUN -eq 1 ]; then
                    echo "[DRY RUN] Would remove: $file (no match within 2 minutes)"
                else
                    echo "Removing unmatched file: $file"
                    rm "$file"
                fi
            fi
        fi
    fi
done

rm archive/temp_nextflow_log

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Pipeline info cleanup simulation complete!"
else
    echo "Pipeline info cleanup complete!"
fi