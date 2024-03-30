
# Documentation for simulating PCR cycles

## Simulation Workflow

This workflow utilizes wf-SLSim to perform sequencing data simulation. It's designed to run in an environment with Nextflow and Docker installed.

### Prerequisites

- Nextflow
- Docker

### Execution

The following script executes the simulation workflow for different PCR cycles:

```bash
#!/bin/bash

# Loop from 1 to 5 for PCR simulations
for i in {1..5}
do
   echo "Running for PCR_$i..."
   nextflow run wf-SLSim/main.nf --projectName PCR_$i \
                                 --outdir results/PCR_$i \
                                 --amp $i \
                                 --barcodes data/PCR/filtered.bc.${i}PCR.csv
done

echo "All runs completed."
```

To execute this script:

1. Save it to a file, e.g., `run_simulation.sh`.
2. Make it executable with `chmod +x run_simulation.sh`.
3. Run the script with `./run_simulation.sh`.

### Output

Results will be stored in separate folders under `results/PCR_x`, where `x` is the PCR simulation cycle.

## Sicelore Workflow

This script uses Docker to run Sicelore. It's designed to run in an environment with Docker installed.

### Prerequisites

- Docker
- Access to fastq sequencing data
- Sicelore Docker image: "genomicpariscentre/sicelore:2.1"

### Execution

The following script executes Sicelore for different PCR cycles:

```bash
#!/bin/bash

# Loop from 1 to 5 for Sicelore analyses
for i in {1..5}
do
   docker run -i --rm \
       -v /export/home1/ScNaUmi-seq_B2022/:/user/local/ScNaUmi-seq_B2022/ \
       -w /user/local/ScNaUmi-seq_B2022/sicelore-2.1/ \
       genomicpariscentre/sicelore:2.1 bash \
       -c "nextflow run sicelore-nf/main.nf \
           --fastqdir /user/local/ScNaUmi-seq_B2022/wf-SLSim/results/PCR_$i/simulated.fastq \
           --PREFIX PCR_$i \
           --project PCR_$i \
           --outdir sicelore-2.1/output/PCR/"
done

echo "All runs completed."
```

To execute this script:

1. Save it to a file, e.g., `run_sicelore.sh`.
2. Make it executable with `chmod +x run_sicelore.sh`.
3. Run the script with `./run_sicelore.sh`.

### Output

Results will be stored in `sicelore-2.1/output/PCR/`, with a sub-folder for each PCR cycle.

---

## General Notes

- Ensure all paths and mounted volumes are correct for your specific environment.
- For long or important workflows, consider using `screen` or `tmux` for a persistent session, or `nohup` for running in the background.

This documentation provides a basic overview to get started with these workflows. Adapt and extend according to the specific needs.
