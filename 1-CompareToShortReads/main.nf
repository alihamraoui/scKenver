
params.samplesheet = './Matrix_samplesheet.csv'

params.outdir = "${projectDir}/output"

process MatrixProcessing {

    input:
    tuple path(scriptDir), path(shortReads), path(longReads), val(dataName), val(dataType)

    output:
    tuple path("${dataName}.rds"), val(dataType), emit: rds
    path 'MatrixProcessing.html', emit: html


    script:
    """
    cp -L $scriptDir/bin/MatrixProcessing.Rmd MatrixProcessing.Rmd

    Rscript -e "rmarkdown::render(
                    'MatrixProcessing.Rmd',
                    params = list(
                        data_type = '${dataType}',
                        data_name = '${dataName}',
                        sr_data = '${shortReads}', 
                        lr_data = '${longReads}')
                    )"
    """
}

process CompareToShortReads {
    input:
    tuple path(scriptDir), path(matrixRds), val(dataType)

    output:
    tuple file('*_QC.per.cell.Rdata'), file('*_total_UMI.Rdata')

    script:
    """
    cp -L $scriptDir/bin/CompareToShortReads/1-processing.Rmd processing.Rmd
    cp -L $scriptDir/bin/imports/common_imports.R common_imports.R
    cp -L $scriptDir/bin/imports/QC.R QC.R

    Rscript -e "rmarkdown::render(
                    'processing.Rmd',
                    params = list(
                        data = '${matrixRds}')
                    )"
    """
}

process VisualizeQC {
    publishDir "${params.outdir}/1-CompareToShortReads/CompareToShortReads", mode: 'copy'
    //cache false

    input:
    path scriptDir
    path files

    output:
    path "QC.pdf"
    path "scatters.pdf"
    path "visualization.html"
    path "correlation_to_SR_metrics.csv", emit: corr_to_SR_metrics 


    script:
    """
    cp -L $scriptDir/bin/CompareToShortReads/2-visualization.Rmd visualization.Rmd
    cp -L $scriptDir/bin/imports/common_imports.R common_imports.R
    cp -L $scriptDir/bin/imports/QC.R QC.R

    Rscript -e "rmarkdown::render(
                    'visualization.Rmd',
                    params = list(
                        files = '${files}')
                    )"
    """
}

process SingleCellMetrics {
    input:
    tuple path(scriptDir), path(matrixRds), val(dataType)

    output:
    tuple file('*_iLISI.Rdata'), file('*_cLISI.Rdata'), file('*_ARI.Rdata')

    script:
    """
    cp -L $scriptDir/bin/singleCellMetrics/scProcessing.Rmd scProcessing.Rmd
    cp -L $scriptDir/bin/imports/common_imports.R common_imports.R
    cp -L $scriptDir/bin/imports/clustering.R clustering.R
    cp -L $scriptDir/bin/imports/cell_annot_custom.R cell_annot_custom.R
    cp -L $scriptDir/bin/imports/clusters_annot.R clusters_annot.R
    cp -L $scriptDir/bin/imports/cell_markers.rda cell_markers.rda
    cp -L $scriptDir/bin/imports/color_markers.rda color_markers.rda

    Rscript -e "rmarkdown::render(
                    'scProcessing.Rmd',
                    params = list(
                        data = '${matrixRds}',
                        nfeatures = '${params.nfeatures}',
                        min_cells = '${params.min_cells}',
                        min_features = '${params.min_features}',
                        reduction = '${params.reduction}',
                        dim = '${params.dim}')
                    )"
    """
}

process ScMetricsVisualisation {
    publishDir "${params.outdir}/singleCellMetrics", mode: 'copy'
    //cache false

    input:
    path scriptDir
    path files

    output:
    path "quanti_metrics.pdf"
    path "scVisualisation.html"
    path "cellTypeMetrics.csv"


    script:
    """
    cp -L $scriptDir/bin/singleCellMetrics/scVisualisation.Rmd scVisualisation.Rmd
    cp -L $scriptDir/bin/imports/common_imports.R common_imports.R
    cp -L $scriptDir/bin/imports/clustering.R clustering.R
    cp -L $scriptDir/bin/imports/cell_markers.rda cell_markers.rda
    cp -L $scriptDir/bin/imports/color_markers.rda color_markers.rda

    Rscript -e "rmarkdown::render(
                    'scVisualisation.Rmd',
                    params = list(
                        files = '${files}')
                    )"
    """
}

process SpatialMetrics {
    publishDir "${params.outdir}/SpatialMetrics", mode: 'copy'

    input:
    tuple path(scriptDir), path(matrixRds), val(dataType)

    output:
    path 'plt.metrics.pdf'
    path 'SpatialDimPlot.pdf'
    path "spatial.csv"

    script:
    """
    cp -L $scriptDir/bin/spatialMetrics/spatial.Rmd spatial.Rmd
    cp -L $scriptDir/bin/imports/spatial.R spatial.R

    Rscript -e "rmarkdown::render(
                    'spatial.Rmd',
                    params = list(
                        data = '${matrixRds}')
                    )"
    """
}

workflow CompareToShortReadsWF{
    take:
        scriptDir
    main:
        Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true, strip: true)
            .map { row ->
                tuple(
                    file(row.shortReads, type: 'dir'),
                    file(row.longReads, type: 'dir'),
                    row.dataName,
                    row.dataType
                )
            }
            .combine(scriptDir)
            .map { shortReads, longReads, dataName, dataType, scriptDir ->
                tuple(scriptDir, shortReads, longReads, dataName, dataType)
            }
            .set { samplesheetTuple }

        MatrixProcessing(samplesheetTuple)

        MatrixProcessing.out.rds
            .combine(scriptDir)
            .map { rds, dataType, scriptDir ->
                tuple(scriptDir, rds, dataType)
            }
            .set { processedTuple }

        CompareToShortReads(processedTuple)

        VisualizeQC(scriptDir, CompareToShortReads.out.collect())

        MatrixProcessing.out.rds
            .filter { it[1].toLowerCase() == 'singlecell' }
            .combine(scriptDir)
            .map { rds, dataType, scriptDir ->
                tuple(scriptDir, rds, dataType)
            }
            .set { singleCellRds }

        SingleCellMetrics(singleCellRds)

        ScMetricsVisualisation(scriptDir, SingleCellMetrics.out.collect())

        MatrixProcessing.out.rds
            .filter { it[1].toLowerCase() == 'spatial' }
            .combine(scriptDir)
            .map { rds, dataType, scriptDir ->
                tuple(scriptDir, rds, dataType)
            }
            .set { spatialRds }
        
        SpatialMetrics(spatialRds)

    emit:
        corr_to_SR_metrics = VisualizeQC.out.corr_to_SR_metrics
        //sc_metrics = ScMetricsVisualisation.out.cellTypeMetrics
        //spatial_metrics = SpatialMetrics.out.spatial
}


workflow {
    scriptDir = Channel.fromPath("${projectDir}")
    CompareToShortReadsWF(scriptDir)
}



