# Benchmark Analysis Methodology: 
## Extracting BC, UMI, and Transcript IDs from Sicelore Data and Comparison with Simulated Gold Standard Data

We focused on evaluating the accuracy and efficiency of isoform detection from single-cell long-read sequencing data. Leveraging the capabilities of Sicelore, we performed a comparison between the assigned labels and the simulated gold standard to validate the accuracy.

This document details the procedures for extracting Barcodes (BC), Unique Molecular Identifiers (UMI), and Transcript IDs (IT) from Sicelore data, followed by their comparison against a gold standard dataset. The steps taken aim to ensure both reproducibility and clarity in our benchmarking process.

### Prerequisites

- Java Runtime Environment (JRE) or Java Development Kit (JDK)
- Sicelore version 2.1

### Data Extraction

We utilized the `jvarkit` toolkit, specifically its `bioalcidaejdk` component, to parse the BAM files produced by Sicelore. This process enabled us to extract the necessary information, including BC, UMI, and IT, directly from each read.
[Visit BioAlcidaeJdk documentation](http://lindenb.github.io/jvarkit/BioAlcidaeJdk.html)
#### Execution Command

```bash
PromethION java -jar /path/to/jvarkit/bioalcidaejdk.jar \
                -e 'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("BC")+"\t"+R.getAttribute("U8")+"\t"+R.getAttribute("IT")));' \
                /path/to/sicelore/03.umis/passedParsed.bam | sed 's/_/\t/5' | cut -f1,3,4 > read_bc_umi_it_sicelore.txt
```

*Please adjust `/path/to/jvarkit/` and `/path/to/sicelore/` to your local directories.*

This command extracts the read names, BC, UMI, and IT for each entry in the BAM file. The output is processed further with `sed` and `cut` to organize the data, facilitating the subsequent comparison phase.

### Data Comparison

Upon extraction, we performed a comprehensive comparison between our dataset and the gold standard. This step was critical for validating the accuracy of the extraction process and ensuring the reliability of the data for downstream analyses.

#### Comparison Procedure

1. **Prepare the Datasets**: Ensure both the extracted dataset and the gold standard are formatted correctly and accessible.
2. **Execute the Comparison**: Utilize a comparison tool or script capable of meticulously analyzing and contrasting the two datasets.
3. **Analyze the Results**: Carefully review the comparison outcomes, noting the accuracy, discrepancies, and any potential anomalies.
