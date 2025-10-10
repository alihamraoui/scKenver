
process SplitGtf {
    input:
    tuple path(gtf), val(tool)

    output:
    tuple path('novel.gtf'), val(tool)

    script:
    """
     python3  $projectDir/bin/splitGtf.py -g $gtf -t $tool -o .
    """
}

process Gffcompare {
    publishDir "${params.outdir}/3-Isoform_discovery/${tool}", mode: 'copy'
    input:
    tuple path(gtf), val(tool), path(truth)

    output:
    path "${tool}_gffcompare.stats"

    script:
    """
    gffcompare -R -r $truth -o ${tool}_gffcompare $gtf
    """
}

process MetricsVisualisation {
    cache false
    publishDir "${params.outdir}/3-Isoform_discovery", mode: 'copy'

    input:
    path scriptDir
    path files

    output:
    path "isoforms_discovery.pdf"
    path "Visualization.html"
    path "Isoforms_discovery_metrics.csv", emit: isoform_metrics


    script:
    """
    cp -L $scriptDir/bin/Visualization.Rmd Visualization.Rmd

    Rscript -e "rmarkdown::render(
                    'Visualization.Rmd',
                    params = list(
                        files = '${files}')
                    )"
    """
}

workflow IsoformDiscovery{
    take:
        scriptDir

    main:
        gtf_novel = Channel.fromPath( params.file ) \
            | splitCsv( header:true, strip:true ) \
            | map { row ->
                tuple( file(row.gtf), row.tool )
            } | SplitGtf

        gtf_tp = Channel.fromPath( params.true_positives ) 


        gtf_novel.map { gtf, tool ->
                tuple( file(gtf), tool, file(params.true_positives) )
            } |  Gffcompare

        MetricsVisualisation(scriptDir, Gffcompare.out.collect())

    emit:
        isoform_metrics = MetricsVisualisation.out.isoform_metrics
    
}

workflow {
    scriptDir = Channel.fromPath("${projectDir}")

    IsoformDiscovery(scriptDir)
}