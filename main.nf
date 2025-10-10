
params.outdir = "${projectDir}/output"

include { CompareToShortReadsWF } from './1-CompareToShortReads/main.nf'
include { UMITranscriptAssignment } from './2-UMI_transcript_assignment/main.nf'
include { IsoformDiscovery } from './3-Isoform_discovery/main.nf'
include { ComparisonWithGroundTruth } from './4-Comparison_with_ground_truth/main.nf'


process MetricsSummary {
    cache false
    publishDir params.outdir, mode:'copy'
    input:
    path metrics_files

    output:
    file 'MetricsSummary.html'
    file 'MetricsSummary.pdf'
    file 'Summary_table.csv'

    script:
    """
    cp -L $projectDir/bin/MetricsSummary.Rmd MetricsSummary.Rmd
    
    Rscript -e "rmarkdown::render('MetricsSummary.Rmd',
                                  params=list(metrics_files='$metrics_files'))"
    """
}


workflow {
    wf1Dir = Channel.fromPath("${projectDir}/1-CompareToShortReads")
    wf2Dir = Channel.fromPath("${projectDir}/2-UMI_transcript_assignment")
    wf3Dir = Channel.fromPath("${projectDir}/3-Isoform_discovery")
    wf4Dir = Channel.fromPath("${projectDir}/4-Comparison_with_ground_truth")

    CompareToShortReadsWF(wf1Dir)
    UMITranscriptAssignment(wf2Dir)
    IsoformDiscovery(wf3Dir)
    ComparisonWithGroundTruth(wf4Dir)

    corr_to_SR_metrics = CompareToShortReadsWF.out.corr_to_SR_metrics.collect()
    //sc_metrics = CompareToShortReadsWF.out.sc_metrics.collect()
    //spatial_metrics = CompareToShortReadsWF.out.spatial_metrics.collect()
    umi_metrics = UMITranscriptAssignment.out.umi_metrics.collect()
    transcript_metrics = UMITranscriptAssignment.out.transcript_metrics.collect()
    iso_metrics = IsoformDiscovery.out.isoform_metrics.collect()
    qm_metrics = ComparisonWithGroundTruth.out.qm_metrics.collect()
    dedup_metrics = ComparisonWithGroundTruth.out.dedup_metrics.collect()
    pb_metrics = ComparisonWithGroundTruth.out.pb_metrics.collect()


    all_metrics = iso_metrics.concat(dedup_metrics)
                    .concat(pb_metrics)
                    .concat(qm_metrics)
                    .concat(umi_metrics)
                    .concat(transcript_metrics)
                    .concat(corr_to_SR_metrics)
                    //.concat(sc_metrics)
                    //.concat(spatial_metrics)

    MetricsSummary(all_metrics.collect()) 
}
