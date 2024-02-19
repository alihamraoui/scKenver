#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    STEP1_readscan()
    STEP1_validbarcodes(STEP1_readscan.out.scancsv)
    STEP2_mapping(STEP1_readscan.out.fastqgz)
    STEP3_umis(STEP2_mapping.out.mappingbam)

    // step 4a (barcoded reads)
    STEP4a_matrix(STEP1_validbarcodes.out.csv, STEP3_umis.out.parsedbam)

    // step 4b (consensus molecules)
    STEP4b_addsequence(STEP3_umis.out.parsedbam, STEP1_readscan.out.fastqgz)
    chrs = STEP4b_getchrs(STEP4b_addsequence.out.parsedbamseq) | splitText | map{it -> it.trim()}
    STEP4b_splitbam(chrs, STEP4b_addsequence.out.parsedbamseq, STEP4b_addsequence.out.parsedbamseqbai) | STEP4b_consensus | STEP4b_concatenate | collectFile | STEP4b_deduplicate | STEP4b_mapping | STEP4b_addtags | STEP4b_addgenes
    STEP4b_matrix(STEP1_validbarcodes.out.csv, STEP4b_addgenes.out.bam)
}

process STEP1_readscan {
    
    output:
    path './passed/ReadScanner.html'        , emit: scanhtml
    path './passed/BarcodesAssigned.tsv'    , emit: scancsv
    path 'fastq_pass.fastq.gz'              , emit: fastqgz

    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
    mkdir ./passed
    $params.java -jar $params.javaXmx $params.nanopore scanfastq -d $params.fastqdir -o ./passed --ncpu $params.max_cpus --bcEditDistance 1 --compress
    find ./passed/passed/ -type f -name '*' | xargs pigz -dc |  pigz > fastq_pass.fastq.gz
    """
}

process STEP1_validbarcodes {
   
    input:
    path(barcodeassigned)

    output:
    path 'BarcodesValidated.csv'    , emit: csv
    
    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
    $params.java -jar $params.javaXmx $params.sicelore SelectValidCellBarcode -I $barcodeassigned -O BarcodesValidated.csv -MINUMI $params.MINUMI -ED0ED1RATIO $params.ED0ED1RATIO
    """
}

process STEP2_mapping {
   
    input:
    path(fastqgz)

    output:
    path 'passed.bam'	, emit: mappingbam
    path 'passed.bam.bai'	, emit: mappingbai
    
    publishDir "${params.outdir}/${params.mappingdir}", mode: 'symlink'
    
    """
    $params.minimap2 -ax splice -uf --sam-hit-only -t $params.max_cpus --junc-bed $params.juncbed $params.minimapfasta $fastqgz | $params.samtools view -bS -@ $params.max_cpus - | $params.samtools sort -m 2G -@ $params.max_cpus -o passed.bam -&& $params.samtools index passed.bam
    """
}

process STEP3_umis {
    
    input:
    path(mappingbam)
 
    output:
    path 'passedParsed.bam'                 , emit: parsedbam
    path 'passedParsed.bai'                 , emit: parsedbai
    path 'passedParsed.bam.genecounts.tsv'  , emit: genecounts
    path 'passedParsed.bam.html'            , emit: umireport
    path 'passedParsed.bam.UMIdepths.tsv'	 , emit: umidepth
    
    publishDir "${params.outdir}/${params.umisdir}", mode: 'copy'
    
    """
    $params.java -jar $params.javaXmx -XX:ActiveProcessorCount=$params.max_cpus $params.nanopore assignumis --inFileNanopore $mappingbam -o passedParsed.bam --annotationFile $params.refflat
    """
}

