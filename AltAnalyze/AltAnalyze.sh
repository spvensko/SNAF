#!/bin/bash

# Species configuration: set SPECIES_PREFIX=Hs (human) or Mm (mouse)
# Set SPECIES_PREFIX and ENSMART_VERSION env vars before calling this script
# For human: SPECIES_PREFIX=Hs, ENSMART_VERSION=EnsMart91 (default)
# For mouse: SPECIES_PREFIX=Mm, ENSMART_VERSION=EnsMart31
SPECIES_PREFIX="${SPECIES_PREFIX:-Hs}"
ENSMART_VERSION="${ENSMART_VERSION:-EnsMart91}"

# Capture work directory before any cd commands
WORK_DIR="$(cd "$(dirname "$0")" && pwd)"

# process the command-line arguments
echo "Current folder is ${WORK_DIR}"
echo "Species prefix: ${SPECIES_PREFIX}, EnsMart version: ${ENSMART_VERSION}"
mode=$1

if [ "$mode" == "bam_to_bed" ]; then
    bam_file=$2
    echo "Running bam to bed workflow, bam file is ${bam_file}"
elif [ "$mode" == "bed_to_junction" ]; then
    bed_folder=$2
    echo "Running bed to junction workflow, bed folder is ${bed_folder}"
elif [ "$mode" == "identify" ]; then
    bam_folder=$2
    cores=$3
    echo "Identify splicing junction, bam folder is ${bam_folder}, using ${cores} cores"
elif [ "$mode" == "DE" ]; then
    output_folder=$2
    group_file=$3
    echo "Identify differentially expressed genes, AltAnalyze output folder is ${output_folder}, group file is ${group_file}"
elif [ "$mode" == "GO" ]; then
    gene_list_file=$2
    echo "Gene enrichment analysis using GO-Elite, gene list file is ${gene_list_file}"
elif [ "$mode" == 'DAS' ]; then
    output_folder=$2
    group_file=$3
    echo "Identify differentailly spliced event, AltAnalyze output folder is ${output_folder}, group file is ${group_file}"
else
    echo "Invalid mode specified"
    exit 1
fi



# bam to bed
if [ "$mode" == "bam_to_bed" ]; then
    function run_BAMtoBED() {

    echo "start to get junction bed"
    python ${WORK_DIR}/altanalyze/import_scripts/BAMtoJunctionBED.py --i $1 \
        --species ${SPECIES_PREFIX} --r ${WORK_DIR}/altanalyze/AltDatabase/${ENSMART_VERSION}/ensembl/${SPECIES_PREFIX}/${SPECIES_PREFIX}_Ensembl_exon.txt

    echo "start to get exon bed"
    python ${WORK_DIR}/altanalyze/import_scripts/BAMtoExonBED.py --i $1  \
        --r ${WORK_DIR}/altanalyze/AltDatabase/${ENSMART_VERSION}/ensembl/${SPECIES_PREFIX}/${SPECIES_PREFIX}.bed --s ${SPECIES_PREFIX}

    return 0
    }

    run_BAMtoBED ${bam_file}

# bed to junction
elif [ "$mode" == "bed_to_junction" ]; then
    # step2: multipath-psi
    task="original"

    ### build necessary folder structure
    mkdir ${WORK_DIR}/altanalyze_output
    mkdir ${WORK_DIR}/altanalyze_output/ExpressionInput

    ### build group file
    touch ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt
    cd ${bed_folder}
    count=0
    for file in *__junction.bed; do
        stream=$(echo $file | sed 's/__junction.bed/.bed/g')
        if [ $(($count%2)) == 0 ]; then
            stream+='\t1\texp'
        else
            stream+='\t2\tctl'
        fi
        echo -e $stream >> ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt
        ((count+=1))
    done
    cd ${WORK_DIR}

    ### build comp file
    touch ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt
    echo -e '1\t2' > ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt

    ### run multipath-psi
    python ${WORK_DIR}/altanalyze/AltAnalyze.py --species ${SPECIES_PREFIX} --platform RNASeq --version ${ENSMART_VERSION} \
        --bedDir ${bed_folder} \
        --output ${WORK_DIR}/altanalyze_output \
        --groupdir ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt \
        --compdir ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
        --runGOElite no

    # step3: process count matrix to only contain PSI junctions
    echo "prune the raw junction count matrix"
    SPECIES_PREFIX="${SPECIES_PREFIX}" python "${WORK_DIR}/prune.py"


# identify
elif [ "$mode" == "identify" ]; then
    # step1: bam to bed
    function run_BAMtoBED() {

    echo "start to get junction bed"
    python ${WORK_DIR}/altanalyze/import_scripts/BAMtoJunctionBED.py --i ${g_bam_folder}/$1 \
        --species ${SPECIES_PREFIX} --r ${WORK_DIR}/altanalyze/AltDatabase/${ENSMART_VERSION}/ensembl/${SPECIES_PREFIX}/${SPECIES_PREFIX}_Ensembl_exon.txt

    echo "start to get exon bed"
    python ${WORK_DIR}/altanalyze/import_scripts/BAMtoExonBED.py --i ${g_bam_folder}/$1  \
        --r ${WORK_DIR}/altanalyze/AltDatabase/${ENSMART_VERSION}/ensembl/${SPECIES_PREFIX}/${SPECIES_PREFIX}.bed --s ${SPECIES_PREFIX}

    return 0
    }

    ### collect for bam file name for parallelization
    cd ${WORK_DIR}/${bam_folder}
    for file in *.bam; do echo $file; done > ${WORK_DIR}/samples.txt
    cd ${WORK_DIR}

    ### start to run
    export -f run_BAMtoBED
    export SPECIES_PREFIX
    export ENSMART_VERSION
    export TMPDIR=/tmp
    export g_bam_folder=${WORK_DIR}/${bam_folder}
    cat ${WORK_DIR}/samples.txt | parallel -P ${cores} run_BAMtoBED {}

    ### move bed files to bed folder
    mkdir ${WORK_DIR}/bed
    cd ${WORK_DIR}/${bam_folder}
    for file in *.bed; do mv $file ${WORK_DIR}/bed; done
    cd ${WORK_DIR}

    # step2: multipath-psi
    task="original"

    ### build necessary folder structure
    mkdir ${WORK_DIR}/altanalyze_output
    mkdir ${WORK_DIR}/altanalyze_output/ExpressionInput

    ### build group file
    touch ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt
    cd ${WORK_DIR}/bed
    count=0
    for file in *__junction.bed; do
        stream=$(echo $file | sed 's/__junction.bed/.bed/g')
        if [ $(($count%2)) == 0 ]; then
            stream+='\t1\texp'
        else
            stream+='\t2\tctl'
        fi
        echo -e $stream >> ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt
        ((count+=1))
    done
    cd ${WORK_DIR}

    ### build comp file
    touch ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt
    echo -e '1\t2' > ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt

    ### run multipath-psi
    python ${WORK_DIR}/altanalyze/AltAnalyze.py --species ${SPECIES_PREFIX} --platform RNASeq --version ${ENSMART_VERSION} \
        --bedDir ${WORK_DIR}/bed \
        --output ${WORK_DIR}/altanalyze_output \
        --groupdir ${WORK_DIR}/altanalyze_output/ExpressionInput/groups.${task}.txt \
        --compdir ${WORK_DIR}/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
        --runGOElite no

    # step3: process count matrix to only contain PSI junctions
    echo "prune the raw junction count matrix"
    SPECIES_PREFIX="${SPECIES_PREFIX}" python "${WORK_DIR}/prune.py"

# DE
elif [ "$mode" == "DE" ]; then
    python ${WORK_DIR}/altanalyze/stats_scripts/metaDataAnalysis.py --p RNASeq --s ${SPECIES_PREFIX} --adjp yes --pval 1 --f 1 \
           --i ${output_folder}/ExpressionInput/exp.original-steady-state.txt \
           --m ${group_file}

# GO
elif [ "$mode" == "GO" ]; then
    # BioMarkers
    mkdir ${WORK_DIR}/GO_Elite_result_BioMarkers
    python ${WORK_DIR}/altanalyze/GO_Elite.py --species ${SPECIES_PREFIX} --mod Ensembl --pval 0.05 --num 3 \
        --input ${gene_list_file} \
        --output ${WORK_DIR}/GO_Elite_result_BioMarkers --dataToAnalyze BioMarkers

    # GO
    mkdir ${WORK_DIR}/GO_Elite_result_GeneOntology
    python ${WORK_DIR}/altanalyze/GO_Elite.py --species ${SPECIES_PREFIX} --mod Ensembl --pval 0.05 --num 3 \
        --input ${gene_list_file} \
        --output ${WORK_DIR}/GO_Elite_result_GeneOntology --dataToAnalyze GeneOntology


# DAS
elif [ "$mode" == 'DAS' ]; then
    python ${WORK_DIR}/altanalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
        --i ${output_folder}/AltResults/AlternativeOutput/${SPECIES_PREFIX}_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
        --m ${group_file}

fi