
params.outdir = "${projectDir}/output"

process TrnsQuantification {
    publishDir params.outdir, mode:'copy'
    input:
    path qm_dir

    output:
    file 'QuantificationMetrics.html'
    file 'QuantificationMetrics.pdf'

    script:
    """
    cp -L $projectDir/bin/Quantification_metrics.Rmd QuantificationMetrics.Rmd
    cp -L $projectDir/bin/imports.R imports.R
    
    Rscript -e "rmarkdown::render('QuantificationMetrics.Rmd',
                                  params=list(data='$qm_dir'))"
    """
}

process UMIDedup {
    publishDir params.outdir, mode:'copy'
    input:
    path dedup_dir

    output:
    file 'UMIDeduplicationMetrics.html'
    file 'UMIDeduplicationMetrics.pdf'

    script:
    """
    cp -L $projectDir/bin/Deduplication_metrics.Rmd UMIDeduplicationMetrics.Rmd
    cp -L $projectDir/bin/imports.R imports.R
    
    Rscript -e "rmarkdown::render('UMIDeduplicationMetrics.Rmd',
                                  params=list(data='$dedup_dir'))"
    """
}

process PseudoBulk {
    publishDir params.outdir, mode:'copy'
    input:
    path qm_dir
    path cell_annotation

    output:
    file 'PseudoBulkMetrics.html'
    file 'Cell_types.pdf'
    file 'PseudoBulkMetrics.pdf'

    script:
    """
    cp -L $projectDir/bin/PseudoBulk_metrics.Rmd PseudoBulkMetrics.Rmd
    cp -L $projectDir/bin/pseudobulk.R pseudobulk.R
    
    Rscript -e "rmarkdown::render('PseudoBulkMetrics.Rmd',
                                  params=list(data='$qm_dir',
                                              cell_annotation='${cell_annotation}'))"
    """
}


workflow {

    qm_dir = Channel.fromPath(params.Qauntification_matrices)
    dedup_dir = Channel.fromPath(params.dedup_matrices)

    cell_annotation_ch = params.cell_annotation != null  ? file(params.cell_annotation) : 
                                             file("no_annotation", type: "file")

    if (cell_annotation_ch.name == "no_annotation") {
            println "\u001B[31mPlease provide the path to the cell annotation file using the '--cell_annotation' option. This is REQUIRED for Pseudobulk metrics analysis.\u001B[0m"
            System.exit(1)
        }

    TrnsQuantification(qm_dir)
    UMIDedup(dedup_dir)
    PseudoBulk(qm_dir, cell_annotation_ch)
}