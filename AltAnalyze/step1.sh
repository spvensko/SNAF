#!/bin/bash


function run() {
    mkdir ./tmp/$1
    cp ./bam/$1.bam ./tmp/$1
    cd ./tmp/$1
    altanalyze bam_to_bed $1.bam
    cd ../..
}

# pre
mkdir tmp
mkdir bed

cd ./bam
for file in *.bam; do echo $(basename -s .bam $file); done > ../samples.txt 
cd ..

# run
export -f run
cat samples.txt | xargs -n 1 -P 4 -I {} bash -c "run {}"

# post
while read line; do 
    mv ./tmp/${line}/${line}.bam.bai ./bam
    mv ./tmp/${line}/*.bed ./bed
    done < samples.txt
rm -r ./tmp
