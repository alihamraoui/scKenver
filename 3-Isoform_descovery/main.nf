
process splitGtf {
    input:
    tuple path(gtf), val(tool)

    output:
    tuple path('novel.gtf'), val(tool)

    script:
    """
     python3  $projectDir/bin/splitGtf.py -g $gtf -t $tool -o .
    """
}

process gffcompare {
    publishDir "${params.outdir}/${tool}", mode: 'copy'
    input:
    tuple path(gtf), val(tool), path(truth)

    output:
    path "${tool}_gffcompare.stats"

    script:
    """
    gffcompare -R -r $truth -o ${tool}_gffcompare $gtf
    """
}

process metricsVisualisation {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path files

    output:
    path "isoforms_discovery.pdf"
    path "Visualization.html"
    path "isoforms_discovery.csv"


    script:
    """
    cp -L $projectDir/bin/Visualization.Rmd Visualization.Rmd

    Rscript -e "rmarkdown::render(
                    'Visualization.Rmd',
                    params = list(
                        files = '${files}')
                    )"
    """
}

workflow { 
    gtf_novel = Channel.fromPath( params.file ) \
        | splitCsv( header:true, strip:true ) \
        | map { row ->
            tuple( file(row.gtf), row.tool )
        } | splitGtf

    gtf_tp = Channel.fromPath( params.true_positives ) 


    gtf_novel.map { gtf, tool ->
            tuple( file(gtf), tool, file(params.true_positives) )
        } |  gffcompare

    metricsVisualisation(gffcompare.out.collect())
}