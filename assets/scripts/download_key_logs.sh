# given the timestamp of the run, download the key logs from the results directory
# note: this script is specific to my local setup ; needs to be run from the LOCAL machine (not the remote server running the pipeline!)
# usage: ./download_key_logs.sh <timestamp>

# get the timestamp of the run
RUN_TIMESTAMP=$1

# if timestand is not provided, parse the latest run timestamp from the pipeline_info directory
# note : the last run file with be in the format of execution_timeline_2026-03-03_13-23-13.html
# note : the format of the timestamp is YYYY-MM-DD_HH-MM-SS
if [ -z "$RUN_TIMESTAMP" ]; then
    RUN_TIMESTAMP=$(ls -td ${BASE_DIR}/results/pipeline_info/* | head -n 1 | cut -d '_' -f 4-5 | cut -d '.' -f 1)
    echo "No timestamp provided. Extracted the latest run timestamp: $RUN_TIMESTAMP"
fi

# Set the base directory
BASE_DIR="/home/Users/pbk1/Somatem"

# archived: download with scp (requires 2 authentications)
# scp -r owl3:${BASE_DIR}/results/pipeline_info/*_${RUN_TIMESTAMP}.* ${BASE_DIR}/.nextflow.log .

# download the key logs 
# note: this will download the pipeline info with the run timestamp and the nextflow log
ssh owl3 "cd ${BASE_DIR} && tar -czf - results/pipeline_info/*_${RUN_TIMESTAMP}.* .nextflow.log" | tar -xzf -

# message about completion  
echo "Key logs downloaded successfully!"

