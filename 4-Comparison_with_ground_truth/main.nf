
params.outdir = "${projectDir}/output"

process TrnsQuantification {
    publishDir "${params.outdir}/4-Comparison_with_ground_truth", mode:'copy'
    input:
    path scriptDir
    path qm_dir

    output:
    file 'QuantificationMetrics.html'
    file 'QuantificationMetrics.pdf'
    path 'Quantification_metrics.csv', emit: qm_metrics

    script:
    """
    cp -L $scriptDir/bin/Quantification_metrics.Rmd QuantificationMetrics.Rmd
    cp -L $scriptDir/bin/imports.R imports.R
    
    Rscript -e "rmarkdown::render('QuantificationMetrics.Rmd',
                                  params=list(data='$qm_dir'))"
    """
}

process UMIDedup {
    publishDir "${params.outdir}/4-Comparison_with_ground_truth", mode:'copy'
    input:
    path scriptDir
    path dedup_dir

    output:
    file 'UMIDeduplicationMetrics.html'
    file 'UMIDeduplicationMetrics.pdf'
    path 'Deduplication_counts_metrics.csv', emit: dedup_metrics

    script:
    """
    cp -L $scriptDir/bin/Deduplication_metrics.Rmd UMIDeduplicationMetrics.Rmd
    cp -L $scriptDir/bin/imports.R imports.R
    
    Rscript -e "rmarkdown::render('UMIDeduplicationMetrics.Rmd',
                                  params=list(data='$dedup_dir'))"
    """
}

process PseudoBulk {
    publishDir "${params.outdir}/4-Comparison_with_ground_truth", mode:'copy'
    input:
    path scriptDir
    path qm_dir
    path cell_annotation

    output:
    file 'PseudoBulkMetrics.html'
    file 'Cell_types.pdf'
    file 'PseudoBulkMetrics.pdf'
    path 'Pseudobulk_metrics.csv', emit: pb_metrics

    script:
    """
    cp -L $scriptDir/bin/PseudoBulk_metrics.Rmd PseudoBulkMetrics.Rmd
    cp -L $scriptDir/bin/pseudobulk.R pseudobulk.R
    
    Rscript -e "rmarkdown::render('PseudoBulkMetrics.Rmd',
                                  params=list(data='$qm_dir',
                                              cell_annotation='${cell_annotation}'))"
    """
}


workflow ComparisonWithGroundTruth{

    take:
        scriptDir

    main:
    qm_dir = Channel.fromPath(params.qauntification_matrices)
    dedup_dir = Channel.fromPath(params.dedup_matrices)

    cell_annotation_ch = params.cell_annotation != null  ? file(params.cell_annotation) : 
                                             file("no_annotation", type: "file")

    if (cell_annotation_ch.name == "no_annotation") {
            println "\u001B[31mPlease provide the path to the cell annotation file using the '--cell_annotation' option. This is REQUIRED for Pseudobulk metrics analysis.\u001B[0m"
            System.exit(1)
        }

    TrnsQuantification(scriptDir, qm_dir)
    UMIDedup(scriptDir, dedup_dir)
    PseudoBulk(scriptDir, qm_dir, cell_annotation_ch)

    emit:
        qm_metrics = TrnsQuantification.out.qm_metrics
        dedup_metrics = UMIDedup.out.dedup_metrics
        pb_metrics = PseudoBulk.out.pb_metrics
}


workflow {
    scriptDir = Channel.fromPath("${projectDir}")

    ComparisonWithGroundTruth(scriptDir)
}
