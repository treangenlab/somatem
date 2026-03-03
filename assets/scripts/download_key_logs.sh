# given the timestamp of the run, download the key logs from the results directory
# note: this script is specific to my local setup ; needs to be run from the LOCAL machine (not the remote server running the pipeline!)
# usage: ./download_key_logs.sh <timestamp>

# get the timestamp of the run
RUN_TIMESTAMP=$1

# Set the base directory
BASE_DIR="/home/Users/pbk1/Somatem"

# archived: download with scp (requires 2 authentications)
# scp -r owl3:${BASE_DIR}/results/pipeline_info/*_${RUN_TIMESTAMP}.* ${BASE_DIR}/.nextflow.log .

# download the key logs 
# note: this will download the pipeline info with the run timestamp and the nextflow log
ssh owl3 "cd ${BASE_DIR} && tar -czf - results/pipeline_info/*_${RUN_TIMESTAMP}.* .nextflow.log" | tar -xzf -

# message about completion  
echo "Key logs downloaded successfully!"

