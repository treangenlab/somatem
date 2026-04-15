#!/usr/bin/env bash

# Better method/ field standard: test on a local build (do before pushing to bioconda)
# runs the full build.sh in an isolated environment and 
# lets you catch substitution/path bugs before pushing to bioconda CI.
# This saves a lot of iteration time on slow CI queues.

# conda build ~/somatem/recipes/somatem/


# ------------------------------------------------------------------
# Use this quick way to test your changes made to the github repo in the bioconda environment
# Author: Windsurf SWE 1.5 fast


name="somatem"
version="0.7.1"

# Copy your updated executable file and update the path to the pipeline directory
cp ~/somatem/bin/somatem ~/micromamba/envs/somatem/bin/somatem
sed -i "s|@@PIPELINE_DIR@@|/home/Users/pbk1/micromamba/envs/somatem/share/${name}-${version}|g" ~/micromamba/envs/somatem/bin/somatem

# make it executable 
chmod +x ~/micromamba/envs/somatem/bin/somatem

# Logs
echo "Updated somatem executable at ~/micromamba/envs/somatem/bin/somatem"
echo "Updated @@PIPELINE_DIR@@ to /home/Users/pbk1/micromamba/envs/somatem/share/${name}-${version}"


# --------------------------------------------------------
# Then symlink the directory for your other changes
# NOTE: This script runs only once for each version of the pipeline to make the symlink and backup
# Comment this out after running once
# --------------------------------------------------------
# cd ~/micromamba/envs/somatem/share/
# mv ${name}-${version} ${name}-${version}.original_backup
# ln -s ~/somatem ${name}-${version}