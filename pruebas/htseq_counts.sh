#bucle para iterar sobre cada archivo bam
for bam in ./*.bam; do
	sample=$(basename "$bam" .bam)
	htseq-count -f bam -r pos -t exon -a 0 \
	"$bam"  ../../gencode.v37lift37.annotation.gtf \
	> ../counts/hisat2/"${sample}_counts.tsv"
done 

