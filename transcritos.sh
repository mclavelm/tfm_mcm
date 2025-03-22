for bam in *.bam; do
	sample=$(basename "$bam" .bam)
	samtools view -H "$bam" | grep '^@SQ' | cut -f2 | sed 's/SN://' >> transcritos_id.txt
done 
#ordenamos y eliminamos duplicados
sort -u transcritos_id.txt -o transcritos_id.txt
