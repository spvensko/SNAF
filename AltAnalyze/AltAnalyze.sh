#!/bin/bash


# process the command-line arguments
echo "Current folder is "$PWD
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
    python ../altanalyze/import_scripts/BAMtoJunctionBED.py --i $1 \
        --species Hs --r ${PWD}/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

    echo "start to get exon bed"
    python ../altanalyze/import_scripts/BAMtoExonBED.py --i $1  \
        --r ${PWD}/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed --s Hs

    return 0
    }  

    run_BAMtoBED ${bam_file}

# bed to junction
elif [ "$mode" == "bed_to_junction" ]; then
    # step2: multipath-psi
    task="original"

    ### build necessary folder structure
    mkdir ${PWD}/altanalyze_output
    mkdir ${PWD}/altanalyze_output/ExpressionInput

    ### build group file
    touch ${PWD}/altanalyze_output/ExpressionInput/groups.${task}.txt
    cd ${bed_folder}
    count=0
    for file in *__junction.bed; do 
        stream=$(echo $file | sed 's/__junction.bed/.bed/g')
        if [ $(($count%2)) == 0 ]; then
            stream+='\t1\texp'
        else
            stream+='\t2\tctl'
        fi
        echo -e $stream >> ../${PWD}/altanalyze_output/ExpressionInput/groups.${task}.txt
        ((count+=1))
    done
    cd ..

    ### build comp file
    touch ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt
    echo -e '1\t2' > ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt

    ### run multipath-psi
    python AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
        --bedDir ${bed_folder} \
        --output ${PWD}/altanalyze_output \
        --groupdir ${PWD}/altanalyze_output/ExpressionInput/groups.${task}.txt \
        --compdir ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
        --runGOElite no

    # step3: process count matrix to only contain PSI junctions
    echo "prune the raw junction count matrix"
    python /usr/src/app/prune.py 


# identify
elif [ "$mode" == "identify" ]; then
    # step1: bam to bed
    function run_BAMtoBED() {
    
    echo "start to get junction bed"
    python ../altanalyze/import_scripts/BAMtoJunctionBED.py --i $1 \
        --species Hs --r ../altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

    echo "start to get exon bed"
    bam_folder=bam
    python ../altanalyze/import_scripts/BAMtoExonBED.py --i $1  \
        --r ../altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed --s Hs

    return 0
    }

    ### collect for bam file name for parallelization
    cd ${bam_folder}
    for file in *.bam; do echo $file; done > ../samples.txt 

    ### start to run
    export -f run_BAMtoBED
    export TMPDIR=/tmp
    export g_bam_folder=${PWD}/${bam_folder}
    cat ../samples.txt | parallel -P ${cores} run_BAMtoBED {}

    ### move bed files to bed folder
    cd ..
    mkdir bed
    cd ${bam_folder}
    for file in *.bed; do mv $file ../bed; done
    cd ..

    # step2: multipath-psi
    task="original"

    ### build necessary folder structure
    mkdir ${PWD}/altanalyze_output
    mkdir ${PWD}/altanalyze_output/ExpressionInput

    ### build group file
    touch ${PWD}/altanalyze_output/ExpressionInput/groups.${task}.txt
    cd ${PWD}/bed
    count=0
    for file in *__junction.bed; do 
        stream=$(echo $file | sed 's/__junction.bed/.bed/g')
        if [ $(($count%2)) == 0 ]; then
            stream+='\t1\texp'
        else
            stream+='\t2\tctl'
        fi
        echo -e $stream >> ../altanalyze_output/ExpressionInput/groups.${task}.txt
        ((count+=1))
    done
    cd ..

    ### build comp file
    touch ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt
    echo -e '1\t2' > ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt

    ### run multipath-psi
    python ${PWD}/altanalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
        --bedDir ${PWD}/bed \
        --output ${PWD}/altanalyze_output \
        --groupdir ${PWD}/altanalyze_output/ExpressionInput/groups.${task}.txt \
        --compdir ${PWD}/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
        --runGOElite no

    # step3: process count matrix to only contain PSI junctions
    echo "prune the raw junction count matrix"
    python /usr/src/app/prune.py 

# DE 
elif [ "$mode" == "DE" ]; then
    python ${PWD}/altanalyze/stats_scripts/metaDataAnalysis.py --p RNASeq --s Hs --adjp yes --pval 1 --f 1 \
           --i ${output_folder}/ExpressionInput/exp.original-steady-state.txt \
           --m ${group_file}

# GO
elif [ "$mode" == "GO" ]; then
    # BioMarkers
    mkdir GO_Elite_result_BioMarkers
    python ${PWD}/altanalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
        --input ${gene_list_file} \
        --output GO_Elite_result_BioMarkers --dataToAnalyze BioMarkers

    # GO
    mkdir GO_Elite_result_GeneOntology
    python ${PWD}/altanalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
        --input ${gene_list_file} \
        --output GO_Elite_result_GeneOntology --dataToAnalyze GeneOntology


# DAS
elif [ "$mode" == 'DAS' ]; then
    python ${PWD}/altanalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
        --i ${output_folder}/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
        --m ${group_file}

fi
