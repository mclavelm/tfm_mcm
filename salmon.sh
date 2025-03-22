#definir directorios
fastq_dir="./"
index_dir="../../salmon_index"
output_dir="../salmon"

#bucle para iterar sobre cada archivo
for fastq in $fastq_dir/*.fastq; do
	sample=$(basename "$fastq" .fastq)
	sample_out="$output_dir/$sample"
	mkdir -p "$sample_out"

	echo "Salmon en muestra $sample"
	salmon quant -i $index_dir -l A -r $fastq --validateMappings -o "$sample_out"
done 
