#Bucle para iterar sobre todos los archivos fastq
for file in trimmed/*.fastq; do
	base=$(basename "$file" _trimmed.fastq)

	#Ejecutar hisat2
	hisat2 -p 1 -x ../hg19_index -U "$file" -S "mapped_out/${base}.sam" \
	2> "mapped_out/${base}_hisat2.log"

	#Convertir SAM a BAM, ordenar e indexar
	samtools sort -m 1G -o "mapped_out/${base}_sorted.bam" "mapped_out/${base}.sam"
	samtools index "mapped_out/${base}_sorted.bam" 

	#Ir eliminando archivos SAM para liberar espacio
	rm "mapped_out/${base}.sam"
done

