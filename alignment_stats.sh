#Definimos los directorios
hisat2_dir="mapped_out_trans"
tmap_dir="bam_iontorrent"
output_file="comparacion_alineamiento.csv"

#Cabecera
echo "sample,total_reads,alignment_rate,mapped_reads,mismatch_rate,aligner" > "$output_file"

#Función para extraer las estadísticas
extract_stats(){
	BAM="$1"
	ALINEADOR="$2"
	SAMPLE=$(basename "$BAM")

	#Total reads
	TOTAL_READS=$(samtools view -c "$BAM")

	#Alignment rate
	ALIGNMENT_RATE=$(samtools flagstat "$BAM" | awk '/mapped \(/ {print $5}' | tr -d '()%')

	#Mapped reads
	MAPPED=$(samtools flagstat "$BAM" | awk '/mapped \(/ {print $1}')

	#Mismatch rate
	MISMATCH=$(samtools stats "$BAM" | awk '/mismatch/ {print $3}' | head -n 1)

	#Añadir resultados a la tabla
	echo "$SAMPLE,$TOTAL_READS,$ALIGNMENT_RATE,$MAPPED,$MISMATCH,$ALINEADOR" >> $output_file
}

#Procesar los BAM de HISAT2
for BAM in $hisat2_dir/*.bam; do
	extract_stats "$BAM" "HISAT2"
done

#Procesar los BAM de TMAP
for BAM in $tmap_dir/*.bam; do
	extract_stats "$BAM" "TMAP"
done 


