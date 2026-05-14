#!/bin/bash

# Script to locate execution_trace_<timestamp>.txt file from the run id within the nextflow log file
# Usage: ./locate_trace_from_nf_log.sh <run_id>

nextflow_log_entries_file="archive/nextflow_log-temp.txt"

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <run_id>"
    exit 1
fi

# Get the run_id
run_id=$1

# locate the timestamp of the run_id within the nextflow log file
timestamp=$(grep "${run_id}" ${nextflow_log_entries_file} | grep "executor" | head -1 | awk '{print $1}')

if [ -z "$timestamp" ]; then
    echo "Timestamp not found for run_id: ${run_id}"
    exit 1
fi

# Find the execution_trace_<timestamp>.txt file within the `results/pipeline_info` directory that is 1 sec before the timestamp
trace_file=$(find . -name "execution_trace_${timestamp}.txt" -path "*/results/pipeline_info/*" 2>/dev/null | head -1)

if [ -z "$trace_file" ]; then
    echo "Trace file not found for run_id: ${run_id}"
    exit 1
fi

# Print the trace file path and contents
echo "Trace file found: \n${trace_file}"
echo "file contents:"
cat ${trace_file}