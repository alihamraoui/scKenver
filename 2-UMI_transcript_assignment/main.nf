
process DataPrep {
    //publishDir params.outdir, mode:'copy'
    input:
    path scriptDir
    path corrected_data
    path raw_data

    output:
    file 'data.rds'

    script:
    """
    cp -L $scriptDir/bin/data_processing.Rmd data_processing.Rmd
    cp -L $scriptDir/bin/imports.r imports.r

    Rscript -e "rmarkdown::render('data_processing.Rmd',
                                    params=list(output_dir='.', 
                                                size=$params.size,
                                                corrected_data='${corrected_data}', 
                                                raw_data='${raw_data}'))"
    """
}

process TranscriptMetrics {
    publishDir "${params.outdir}/2-UMI_transcript_assignment", mode:'copy'
    input:
    path scriptDir
    file data

    output:
    file 'transcriptMetrics.html'
    file 'transcript_assignement.pdf'
    path 'transcript_assignment_metrics.csv', emit: transcript_metrics

    script:
    """
    cp -L $scriptDir/bin/transcript.Rmd transcriptMetrics.Rmd
    cp -L $scriptDir/bin/imports.r imports.r
    
    Rscript -e "rmarkdown::render('transcriptMetrics.Rmd',
                                  params=list(data='$data',
                                              size=$params.size))"
    """
}

process UMIMetrics {
    publishDir "${params.outdir}/2-UMI_transcript_assignment", mode:'copy'
    input:
    path scriptDir
    file data

    output:
    file 'umiMetrics.html'
    file 'UMI_deduplication.pdf'
    file 'UMI_deduplication.pdf'
    path 'UMI_correction_metrics.csv', emit: umi_metrics

    script:
    """
    cp -L $scriptDir/bin/umi.Rmd umiMetrics.Rmd
    cp -L $scriptDir/bin/imports.r imports.r
    
    Rscript -e "rmarkdown::render('umiMetrics.Rmd',
                                  params=list(data='$data',
                                              size=$params.size))"
    """
}

workflow UMITranscriptAssignment{
    take:
        scriptDir
    main:
    corrected_data = Channel.fromPath(params.data_corrected)
    raw_data =  Channel.fromPath(params.data_raw)

    DataPrep(scriptDir, corrected_data, raw_data)
    TranscriptMetrics(scriptDir, DataPrep.out)
    UMIMetrics(scriptDir, DataPrep.out)

    emit:
        umi_metrics = UMIMetrics.out.umi_metrics
        transcript_metrics = TranscriptMetrics.out.transcript_metrics
}

workflow {
    scriptDir = Channel.fromPath("${projectDir}")

    UMITranscriptAssignment(scriptDir)
}
