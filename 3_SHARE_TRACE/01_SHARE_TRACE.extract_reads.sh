#!/bin/bash
R1path="" ## path to fastq for R1 of barcode library

myPATH="../tools/"
matchSeq1="TAGACAT"	# vector sequence to match for extracting reads
numMM="1"		# mismatch to allow in vector extracttion
ncores="4"		# parallel cores for Rscript
twoRead="T"		# two fastq format ("T" or "F")

prefix=$(echo $R1path | sed -e 's/.R1.fastq.gz//')
output=$(echo $R1path | sed -e's/.*\///' | sed -e 's/.R1.fastq.gz//')
date

echo "Extracting vector containing reads"
if [ $twoRead = "T" ]; then
        python3 -u $myPATH/trace_share_01.match_vector.v1.1.py \
                -R1 $prefix.R1.fastq.gz \
                -R2 $prefix.R2.fastq.gz \
                -v $matchSeq1 \
                -m $numMM \
                --out $output
else
	python3 -u $myPATH/trace_share_01.match_vector.v1.1.py \
		-R1 $prefix.R1.fastq.gz \
		-R2 $prefix.I1.fastq.gz \
		-R3 $prefix.I2.fastq.gz \
		-R4 $prefix.R2.fastq.gz \
		-v $matchSeq1 \
		-m $numMM \
		--out $output
fi
date

