
process dataPrep {
    publishDir params.outdir, mode:'copy'
    input:
    path corrected_data
    path raw_data

    output:
    file 'data.rds'

    script:
    """
    cp -L $projectDir/bin/data_processing.Rmd data_processing.Rmd
    cp -L $projectDir/bin/imports.r imports.r

    Rscript -e "rmarkdown::render('data_processing.Rmd',
                                    params=list(output_dir='.', 
                                                size=$params.size,
                                                corrected_data='${corrected_data}', 
                                                raw_data='${raw_data}'))"
    """
}

process transcriptMetrics {
    publishDir params.outdir, mode:'copy'
    input:
    file data

    output:
    file 'transcriptMetrics.html'
    file 'transcript_assignement.pdf'
    file 'transcript_assignement.csv'

    script:
    """
    cp -L $projectDir/bin/transcript.Rmd transcriptMetrics.Rmd
    cp -L $projectDir/bin/imports.r imports.r
    
    Rscript -e "rmarkdown::render('transcriptMetrics.Rmd',
                                  params=list(data='$data',
                                              size=$params.size))"
    """
}

process UMIMetrics {
    publishDir params.outdir, mode:'copy'
    input:
    file data

    output:
    file 'umiMetrics.html'
    file 'UMI_deduplication.pdf'
    file 'UMI_deduplication.pdf'

    script:
    """
    cp -L $projectDir/bin/umi.Rmd umiMetrics.Rmd
    cp -L $projectDir/bin/imports.r imports.r
    
    Rscript -e "rmarkdown::render('umiMetrics.Rmd',
                                  params=list(data='$data',
                                              size=$params.size))"
    """
}

workflow { 
    corrected_data = Channel.fromPath(params.data_corrected)
    raw_data =  Channel.fromPath(params.data_raw)

    dataPrep(corrected_data, raw_data)
    transcriptMetrics(dataPrep.out)
    UMIMetrics(dataPrep.out)
}
