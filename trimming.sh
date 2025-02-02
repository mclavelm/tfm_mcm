#Bucle para procesar cada archivo fastq
for file in *.fastq; do
  #Extraer el nombre base del archivo
  base=$(basename "$file" .fastq)

  #Ejecutar fastp para el archivo 
  fastp -i "$file" \
        -o "../trimmed/${base}_trimmed.fastq"  \
        -q 20 --cut_tail --length_required 100

done

