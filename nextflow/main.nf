#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

params.outdir = "output/"
params.start_dir = "${launchDir}"


Channel.fromPath( "samplesheet.csv" )
        .splitCsv(header: false)
        .set{ samples }
Channel.fromPath("1_processing.Rmd")
        .set{ analysis_ch }

workflow {
    Processing(samples, analysis_ch, params.start_dir)
}

process Processing { 
    input:
    tuple val(sequencing), val(protocole) 
    path(analysis)
    val(data)

    output:
    path '*', emit: processed_data

    publishDir "${params.outdir}", mode: 'copy'
    
    """
    mkdir -p html/
    Rscript -e "rmarkdown::render(input = '${analysis}', 
                                      output_file = 'html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}',
                                      benchmark_dir = '${data}'))"
                                      
    """
}

