#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

params.outdir = "output/"


Channel.fromPath( "test/samplesheet.csv" )
        .splitCsv(header: false)
        .set{ samples }

workflow {
    Processing(samples)
}

process Processing { 
    input:
    tuple val(sequencing), val(protocole) 

    output:
    path 'output/data'        , emit: processed_data

    publishDir "${params.outdir}", mode: 'copy'
    
    """
    Rscript -e "rmarkdown::render(input = '1-processing.Rmd', 
                                      output_file = 'html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}'))"
    """
}