!/bin/bash

# Check if main.nf exists
if [ ! -f "main.nf" ]; then
    echo "ERROR: cannot found main.nf." >&2
    exit 1
fi

# Check if wget is installed
if [ ! -x "$(command -v wget)" ]; then
    echo "ERROR: wget is not installed/not in the PATH." >&2
    exit 1
fi

# Download dataset
if [ ! -f "test_data.zip" ]; then
    wget https://zenodo.org/records/17432963/files/test_data.zip
fi

# Unzip dataset
if [ ! -d "test_data" ]; then
    unzip test_data.zip
fi

# Check if nextflow is in the path
if [ ! -x "$(command -v nextflow)" ]; then
    echo "ERROR: nextflow is not installed/not in the PATH." >&2
    exit 1
fi

echo "* Execute basic workflow"
fasle #nextflow run main.nf --outdir results_global   -with-report -with-trace -with-timeline
