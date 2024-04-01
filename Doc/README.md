# Benchmark Analysis Methodology: 
## Extracting BC, UMI, and Transcript IDs from Sicelore Data and Comparison with Simulated Gold Standard Data

We focused on evaluating the accuracy and efficiency of isoform detection from single-cell long-read sequencing data. Leveraging the capabilities of Sicelore, we performed a comparison between the assigned labels and the simulated gold standard to validate the accuracy.

This document details the procedures for extracting Barcodes (BC), Unique Molecular Identifiers (UMI), and Transcript IDs (IT) from Sicelore data, followed by their comparison against a gold standard dataset. The steps taken aim to ensure both reproducibility and clarity in our benchmarking process.

### Prerequisites

- Java Runtime Environment (JRE) or Java Development Kit (JDK)
- Sicelore version 2.1

## Data Extraction and Comparison

### Sicelore

For Sicelore we employ the  `bioalcidaejdk` component of `jvarkit` toolkit, to parse the BAM files produced by Sicelore. This process enabled us to extract the necessary information from `04a.matrices/isobam.bam`, including BC, UMI, and IT, directly from each read.
[Visit BioAlcidaeJdk documentation](http://lindenb.github.io/jvarkit/BioAlcidaeJdk.html)

```bash
PromethION java -jar /path/to/jvarkit/bioalcidaejdk.jar \
                -e 'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("BC")+"\t"+R.getAttribute("U8")+"\t"+R.getAttribute("IT")));' \
                /path/to/sicelore/04a.matrices/isobam.bam | sed 's/_/\t/4'| cut -f1,3,4,5 | grep -v undef > read_bc_umi_trns_sicelore.tsv
```

### FLAMES

For FLAMES, data is extracted from the `realign2transcript.bam` file using the following command:

```bash
samtools view realign2transcript.bam | sed 's/#/	/1; s/_/	/1'| awk '{print $3 "	" $1 "	" $2 "	" $5}' > read_bc_umi_trns_flames.tsv
```

This command organizes the extracted data into a TSV file, ready for comparison with the gold standard.

## Comparison with Simulated Gold Standard Data

After data extraction, a thorough comparison with simulated gold standard data is performed for each tool. This crucial step validates the accuracy of the extraction process and ensures the reliability of the data for downstream analyses.

#### Comparison Procedure

1. **Prepare the Datasets**: Ensure both the extracted dataset and the gold standard are formatted correctly and accessible.
2. **Execute the Comparison**: Utilize a comparison tool or script capable of meticulously analyzing and contrasting the two datasets.
3. **Analyze the Results**: Carefully review the comparison outcomes, noting the accuracy, discrepancies, and any potential anomalies.
