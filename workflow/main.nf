#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

params.processing = "output/1_processing"

params.start_dir = "${launchDir}"


Channel.fromPath( "samplesheet.csv" )
        .splitCsv(header: false)
        .map { it.join("-") }
        .set{ protocole_ch }

Channel.fromPath("./Rmarkdowns/1_processing.Rmd")
        .set{ processing_Rmd_ch }

Channel.fromPath("./Rmarkdowns/2_visualization.Rmd")
        .set{ visualization_Rmd_ch }

workflow {
    Processing(protocole_ch, processing_Rmd_ch, params.start_dir)
}

process Processing { 
    input:
    each protocole
    path(processing_Rmd)
    val(data)

    output:
    path "html/${protocole}.html", emit: process_html
    path "data/*.Rdata", emit: processed_data

    publishDir "${params.processing}", mode: 'copy'
    
    """
    mkdir -p html/
    mkdir -p data/
    Rscript -e "rmarkdown::render(input = '${processing_Rmd}', 
                                      output_file = 'html/${protocole}.html', 
                                      params = list(protocole = '${protocole}',
                                                    benchmark_dir = '${data}'))"
                                      
    """
}

