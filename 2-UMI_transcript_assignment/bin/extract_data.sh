# Description: Extracts the barcode and UMI from a bam file
# Usage: bash extract_data.sh -b <bam> -o <output> -t <tool>
# Dependencies: JAVA; bioalcidaejdk.jar
# Date: 2021-07-29

while getopts "b:o:t:p:" opt; do
    case $opt in
        b) bam=$OPTARG;;
        o) output=$OPTARG;;
        t) tool=$OPTARG;;
        p) prefix=$OPTARG;;
        \?) echo "Invalid option: $OPTARG";;
    esac
done

if [ -z $bam ] || [ -z $output ] || [ -z $tool ]; then
    echo "Usage: $0 -b <bam> -o <output> -t <tool>"
    exit 1
fi

if [ ! -f $bam ]; then
    echo "File $bam not found"
    exit 1
fi

if [ $tool != "Sicelore" ] && [ $tool != "Sockeye" ] && [ $tool != "FLAMES" ] && [ $tool != "bambu" ]; then
    echo "Tool $tool not recognized"
    exit 1
fi

if [ ! -d ${output}/corrected ]; then
    mkdir ${output}/corrected
fi

if [ ! -d ${output}/raw ]; then
    mkdir ${output}/raw
fi


mkdir -p ${output}/corrected/$tool

mkdir -p ${output}/raw/$tool


if [ $tool == "Sicelore" ]; then
    java -jar bin/dist/bioalcidaejdk.jar -e  'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("BC")+"\t"+R.getAttribute("U8")+"\t"+R.getAttribute("IT")));' $bam | sed 's/_/\t/'| cut -f1,3,4,5 | grep -v undef > ${output}/corrected/Sicelore/Sicelore_$prefix.tsv
    java -jar bin/dist/bioalcidaejdk.jar -e  'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("U7")));' $bam | sed 's/_/\t/' | cut -f1,3,4 > $output/raw/Sicelore/Sicelore_$prefix.tsv
fi

if [ $tool == "Sockeye" ]; then
    java -jar bin/dist/bioalcidaejdk.jar -e  'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("CB")+"\t"+R.getAttribute("UB")+"\t"+R.getAttribute("TR")));' $bam | sed 's/_/\t/'| cut -f1,3,4,5 | grep -Pv '\t-\s*$' > ${output}/corrected/Sockeye/Sockeye_$prefix.tsv
    java -jar bin/dist/bioalcidaejdk.jar -e  'stream().forEach(R->println(R.getReadName()+"\t"+R.getAttribute("UR")));' $bam | sed 's/_/\t/' | cut -f1,4 > $output/raw/Sockeye/Sockeye_$prefix.tsv
fi

if [ $tool == "FLAMES" ]; then
    refgtf='/export/home1/ScNaUmi-seq_B2022/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    mkdir -p tmp/

    echo "extracting barcode and transcript from $bam"
    samtools view $bam | sed 's/_/\t/; s/#/\t/' | awk -F'\t' '{print $3"\t"$1"\t"$5}' | grep -v ENSG > tmp/trx_FLAMES_$prefix.tsv

    echo "extracting UMI from $bam"
    python3 /export/home1/ScNaUmi-seq_B2022/tmp/FLAMES/FLAMES/python/report_umi.py  -a  $refgtf  -f tmp/FLAMES_$prefix.tsv --outdir ${bam%/*}

    echo "Joining barcode, UMI and transcript, sorting ..."
    cut -f1,3  tmp/FLAMES_$prefix.tsv > tmp/umi_FLAMES_$prefix.tsv

    sort -k1,1 tmp/umi_FLAMES_$prefix.tsv > tmp/sorted_corrected_FLAMES_$prefix.tsv
    sort -k1,1 tmp/trx_FLAMES_$prefix.tsv > tmp/sorted_trx_FLAMES_$prefix.tsv

    echo "Joining barcode, UMI and transcript, Joining ..."
    join -t $'\t' tmp/sorted_corrected_FLAMES_$prefix.tsv tmp/sorted_trx_FLAMES_$prefix.tsv | awk -F'\t' '{ print $1"\t"$3"\t"$2"\t"$4}' > $output/corrected/FLAMES/FLAMES_$prefix.tsv

    echo "Save raw UMI"
    cut -f1,2  tmp/FLAMES_$prefix.tsv > $output/raw/FLAMES/FLAMES_$prefix.tsv

    rm -rf tmp/
fi

if [ $tool == "bambu" ]; then
    samtools view $bam | cut -f1 | cut -d '_' -f1,2 | sed 's/#/\t/g;s/_/\t/g' | awk '{print $3"\t"$1"\t"$2"\tENSTUNKNOWN"}' > $output/corrected/bambu/bambu_$prefix.tsv
    cut -f1,3 $output/corrected/bambu/bambu_$prefix.tsv > $output/raw/bambu/bambu_$prefix.tsv
fi