#!/bin/bash

prefix="mouse_colitis_tissue"
fragFile="frags.list" # list of fragment files to use
Genome="hg38"
large="T"

cores=4
myPATH="../tools/"

nfrags=$(wc -l $fragFile | awk '{print $1}')

if [ $nfrags -gt 1 ]; then
	echo "More than one fragment file detected -- merging..."
	if [ $large == "T" ]; then
		echo "Large files -- processing in stepwise fashion"
		zcat $(cat $fragFile) | pigz --fast -p $cores > $prefix.fragments.bed.gz
		zcat $prefix.fragments.bed.gz | sort -k1,1 -k2,2n | bgzip > $prefix.fragments.tsv.gz
		rm $prefix.fragments.bed.gz
	else
		zcat $(cat $fragFile) | bedtools sort -i stdin | bgzip > $prefix.fragments.tsv.gz
	fi
	echo "Indexing"
	tabix -p bed $prefix.fragments.tsv.gz
	bedFile="$prefix.fragments.tsv.gz"
else
	bedFile=$(cat $fragFile)
fi


echo "Calling peaks"
if [ $Genome == 'hg19' ] || [ $Genome == 'hg38' ]; then
        G="hs"
elif [ $Genome == 'mm10' ]; then
        G="mm"
fi

macs2 callpeak -t $bedFile -f BED -g $G -n $prefix.macs2 \
        --nomodel --nolambda --keep-dup all --call-summits \
        --outdir ./

echo "Cleaning summit file" ## for creating a peak set that is of the same width and centered at peak summits
Rscript $myPATH/cleanSummits.v1.1.R $prefix $Genome

echo "Done!"
