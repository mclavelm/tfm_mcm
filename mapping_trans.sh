#Bucle para iterar sobre todos los archivos fastq
for file in trimmed/*.fastq; do
	base=$(basename "$file" _trimmed.fastq)

	#Ejecutar hisat2
	hisat2 -p 1 -x ../transcriptome_index -U "$file" -S "mapped_out_trans/${base}_trans.sam" \
	2> "mapped_out_trans/${base}_hisat2.log"

	#Convertir SAM a BAM, ordenar e indexar
	samtools sort -m 1G -o "mapped_out_trans/${base}_sorted_trans.bam" "mapped_out_trans/${base}_trans.sam"
	samtools index "mapped_out_trans/${base}_sorted_trans.bam"

	#Ir eliminando archivos SAM para liberar espacio
	rm "mapped_out_trans/${base}_trans.sam"

done 
