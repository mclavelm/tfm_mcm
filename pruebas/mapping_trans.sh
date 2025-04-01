#Bucle para iterar sobre todos los archivos fastq
for file in trimmed/*.fastq; do
	base=$(basename "$file" _trimmed.fastq)

	#Ejecutar hisat2
	hisat2 -p 1 -x ../transcriptome_index -U "$file" -S "mapped_out_trans_max1/${base}_trans_max.sam" -k 1 \
	2> "mapped_out_trans_max1/${base}_hisat2_max.log"

	#Convertir SAM a BAM, ordenar e indexar
	samtools sort -m 1G -o "mapped_out_trans_max1/${base}_sorted_trans_max.bam" "mapped_out_trans_max1/${base}_trans_max.sam"
	samtools index "mapped_out_trans_max1/${base}_sorted_trans_max.bam"

	#Ir eliminando archivos SAM para liberar espacio
	rm "mapped_out_trans_max1/${base}_trans_max.sam"

done 
