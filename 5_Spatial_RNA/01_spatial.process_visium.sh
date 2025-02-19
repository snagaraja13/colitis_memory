batch=""
bclPath="" # path to BCLs
flowcell="HJGL2DRX5"

csvPath="" # Lane,Sample,Index. e.g. *,APC01-A1-13-B,SI-TS-D1
sampleFile="" # info for 10x processing: APC01-A1-13-B,A1,H1-MY2FRKC
imageList="" # list of paths to tif images from Visium HD

genome="mm10"

dir="" # raw data diectory

## Count gene expresson
fastqPath="$dir/$batch/outs/fastq_path/$flowcell/"
## reference files from 10x website
if [[ $genome == "mm10" ]]; then
	probeSet="/home/snagaraja/software/spaceranger-3.0.1/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv"
	genomeRef="/mnt/users/snagaraja/ref/spaceranger/refdata-gex-mm10-2020-A/"
elif [[ $genome == "hg38" ]]; then
	probeSet="/home/snagaraja/software/spaceranger-3.0.1/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
	genomeRef="/mnt/users/snagaraja/ref/spaceranger/refdata-gex-GRCh38-2020-A/"
else
	echo "Genome probes and reference files not found"
fi

i=1
while read line; do 
	sample=$(echo $line | sed 's/,.*//g')
	area=$(echo $line | sed 's/^[^,]*,//' | sed 's/,.*//g')
	slide=$(echo $line | sed 's/.*,//g')
	imageFile=$(head -${i} $imageList | tail -1)
	
	echo $sample
	spaceranger count --id=${batch}_${sample} \
		--probe-set=$probeSet \
		--cytaimage=$imageFile \
		--create-bam=true \
		--fastqs=$fastqPath \
		--sample=$sample \
		--transcriptome=$genomeRef \
		--slide=$slide \
		--area=$area
	
	i=$((i+1))
done < $sampleFile

echo "Done"


