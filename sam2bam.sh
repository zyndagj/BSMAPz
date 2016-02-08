#convering SAM to BAM, sort and index BAM
#usage ./sam2bam.sh <infile>
#outputfile will be <infile_stem>.bam and <infile_stem>.bam.bai
#rrbsmap_path=${0%/*}
tmpbam=${1%.*}.tmp.bam
outbam=${1%.*}.bam

echo "Converting SAM to BAM ..."
if [ ! -f $1 ]; then
	echo "$1 does not exist."
	exit 1
fi
samtools/samtools view -bS $1 > $tmpbam
if [ $? -ne 0 ]; then
	echo "SAM2BAM conversion not sucessful."
	echo "$1 remains unchanged."
	rm $tmpbam
	exit 1
fi
echo "Sorting BAM ..."
samtools/samtools sort $tmpbam ${outbam%.*}
if [ $? -ne 0 ]; then
	echo "BAM file sorting not sucessful."
	echo "$outbam is in unsorted BAM format".
	mv $tmpbam $outbam
	exit 1
fi
rm $tmpbam
echo "Indexing BAM ..."
samtools/samtools index $outbam
if [ $? -ne 0 ]; then
	echo "BAM file indexing not sucessful."
	exit 1
fi
exit 0
