#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    Processing()
    Visialization(Processing.out.processed_data)
}

process Processing { 
    
    output:
    path './output/data'        , emit: processed_data

    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
    mkdir -p ./output
    Rscript -e "rmarkdown::render(input = '1-processing.Rmd', 
                                      output_file = 'output/html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}'))"
    """
}

process Visialization { 
    
    output:
    path './output/visualization.html'        , emit: vizhtml


    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
    Rscript -e "rmarkdown::render(input='2-visualization.Rmd', output_file='output/html/visualization.html')"
    """
}
