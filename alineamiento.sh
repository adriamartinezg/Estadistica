#!/bin/bash
# PRJNA701955

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
gzip -d GCF_000002985.6_WBcel235_genomic.fna.gz
mv GCF_000002985.6_WBcel235_genomic.fna celegans.fna

echo "Indexando el genoma"
bowtie2-build celegans.fna celegans


wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/021/SRR13712621/SRR13712621_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/021/SRR13712621/SRR13712621_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712621_1.fastq.gz -2 SRR13712621_2.fastq.gz -S SRR13712621.sam
rm SRR13712621_1.fastq.gz SRR13712621_2.fastq.gz
samtools view -bS SRR13712621.sam > SRR13712621.bam
rm SRR13712621.sam
samtools sort SRR13712621.bam -o SRR13712621.sorted.bam
rm SRR13712621.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/022/SRR13712622/SRR13712622_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/022/SRR13712622/SRR13712622_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712622_1.fastq.gz -2 SRR13712622_2.fastq.gz -S SRR13712622.sam
rm SRR13712622_1.fastq.gz SRR13712622_2.fastq.gz
samtools view -bS SRR13712622.sam > SRR13712622.bam
rm SRR13712622.sam
samtools sort SRR13712622.bam -o SRR13712622.sorted.bam
rm SRR13712622.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/023/SRR13712623/SRR13712623_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/023/SRR13712623/SRR13712623_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712623_1.fastq.gz -2 SRR13712623_2.fastq.gz -S SRR13712623.sam
rm SRR13712623_1.fastq.gz SRR13712623_2.fastq.gz
samtools view -bS SRR13712623.sam > SRR13712623.bam
rm SRR13712623.sam
samtools sort SRR13712623.bam -o SRR13712623.sorted.bam
rm SRR13712623.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/024/SRR13712624/SRR13712624_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/024/SRR13712624/SRR13712624_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712624_1.fastq.gz -2 SRR13712624_2.fastq.gz -S SRR13712624.sam
rm SRR13712624_1.fastq.gz SRR13712624_2.fastq.gz
samtools view -bS SRR13712624.sam > SRR13712624.bam
rm SRR13712624.sam
samtools sort SRR13712624.bam -o SRR13712624.sorted.bam
rm SRR13712624.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/025/SRR13712625/SRR13712625_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/025/SRR13712625/SRR13712625_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712625_1.fastq.gz -2 SRR13712625_2.fastq.gz -S SRR13712625.sam
rm SRR13712625_1.fastq.gz SRR13712625_2.fastq.gz
samtools view -bS SRR13712625.sam > SRR13712625.bam
rm SRR13712625.sam
samtools sort SRR13712625.bam -o SRR13712625.sorted.bam
rm SRR13712625.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/026/SRR13712626/SRR13712626_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/026/SRR13712626/SRR13712626_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712626_1.fastq.gz -2 SRR13712626_2.fastq.gz -S SRR13712626.sam
rm SRR13712626_1.fastq.gz SRR13712626_2.fastq.gz
samtools view -bS SRR13712626.sam > SRR13712626.bam
rm SRR13712626.sam
samtools sort SRR13712626.bam -o SRR13712626.sorted.bam
rm SRR13712626.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/027/SRR13712627/SRR13712627_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/027/SRR13712627/SRR13712627_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712627_1.fastq.gz -2 SRR13712627_2.fastq.gz -S SRR13712627.sam
rm SRR13712627_1.fastq.gz SRR13712627_2.fastq.gz
samtools view -bS SRR13712627.sam > SRR13712627.bam
rm SRR13712627.sam
samtools sort SRR13712627.bam -o SRR13712627.sorted.bam
rm SRR13712627.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/028/SRR13712628/SRR13712628_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/028/SRR13712628/SRR13712628_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712628_1.fastq.gz -2 SRR13712628_2.fastq.gz -S SRR13712628.sam
rm SRR13712628_1.fastq.gz SRR13712628_2.fastq.gz
samtools view -bS SRR13712628.sam > SRR13712628.bam
rm SRR13712628.sam
samtools sort SRR13712628.bam -o SRR13712628.sorted.bam
rm SRR13712628.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/029/SRR13712629/SRR13712629_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/029/SRR13712629/SRR13712629_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712629_1.fastq.gz -2 SRR13712624_9.fastq.gz -S SRR13712629.sam
rm SRR13712629_1.fastq.gz SRR13712629_2.fastq.gz
samtools view -bS SRR13712629.sam > SRR13712629.bam
rm SRR13712629.sam
samtools sort SRR13712629.bam -o SRR13712629.sorted.bam
rm SRR13712629.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/031/SRR13712631/SRR13712631_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/031/SRR13712631/SRR13712631_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712631_1.fastq.gz -2 SRR13712631_2.fastq.gz -S SRR13712631.sam
rm SRR13712631_1.fastq.gz SRR13712631_2.fastq.gz
samtools view -bS SRR13712631.sam > SRR13712631.bam
rm SRR13712631.sam
samtools sort SRR13712631.bam -o SRR13712631.sorted.bam
rm SRR13712631.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/030/SRR13712630/SRR13712630_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/030/SRR13712630/SRR13712630_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712630_1.fastq.gz -2 SRR13712630_2.fastq.gz -S SRR13712630.sam
rm SRR13712630_1.fastq.gz SRR13712630_2.fastq.gz
samtools view -bS SRR13712630.sam > SRR13712630.bam
rm SRR13712630.sam
samtools sort SRR13712630.bam -o SRR13712630.sorted.bam
rm SRR13712630.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/032/SRR13712632/SRR13712632_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/032/SRR13712632/SRR13712632_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712632_1.fastq.gz -2 SRR13712632_2.fastq.gz -S SRR13712632.sam
rm SRR13712632_1.fastq.gz SRR13712632_2.fastq.gz
samtools view -bS SRR13712632.sam > SRR13712632.bam
rm SRR13712632.sam
samtools sort SRR13712632.bam -o SRR13712632.sorted.bam
rm SRR13712632.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/033/SRR13712633/SRR13712633_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/033/SRR13712633/SRR13712633_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712633_1.fastq.gz -2 SRR13712633_2.fastq.gz -S SRR13712633.sam
rm SRR13712633_1.fastq.gz SRR13712633_2.fastq.gz
samtools view -bS SRR13712633.sam > SRR13712633.bam
rm SRR13712633.sam
samtools sort SRR13712633.bam -o SRR13712633.sorted.bam
rm SRR13712633.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/034/SRR13712634/SRR13712634_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/034/SRR13712634/SRR13712634_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712634_1.fastq.gz -2 SRR13712634_2.fastq.gz -S SRR13712634.sam
rm SRR13712634_1.fastq.gz SRR13712634_2.fastq.gz
samtools view -bS SRR13712634.sam > SRR13712634.bam
rm SRR13712634.sam
samtools sort SRR13712634.bam -o SRR13712634.sorted.bam
rm SRR13712634.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/035/SRR13712635/SRR13712635_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/035/SRR13712635/SRR13712635_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712635_1.fastq.gz -2 SRR13712635_2.fastq.gz -S SRR13712635.sam
rm SRR13712635_1.fastq.gz SRR13712635_2.fastq.gz
samtools view -bS SRR13712635.sam > SRR13712635.bam
rm SRR13712635.sam
samtools sort SRR13712635.bam -o SRR13712635.sorted.bam
rm SRR13712635.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/036/SRR13712636/SRR13712636_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/036/SRR13712636/SRR13712636_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712636_1.fastq.gz -2 SRR13712636_2.fastq.gz -S SRR13712636.sam
rm SRR13712636_1.fastq.gz SRR13712636_2.fastq.gz
samtools view -bS SRR13712636.sam > SRR1371236.bam
rm SRR13712636.sam
samtools sort SRR13712636.bam -o SRR13712636.sorted.bam
rm SRR13712636.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/037/SRR13712637/SRR13712637_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/037/SRR13712637/SRR13712637_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712637_1.fastq.gz -2 SRR13712637_2.fastq.gz -S SRR13712637.sam
rm SRR13712637_1.fastq.gz SRR13712637_2.fastq.gz
samtools view -bS SRR13712637.sam > SRR13712637.bam
rm SRR13712637.sam
samtools sort SRR13712637.bam -o SRR13712637.sorted.bam
rm SRR13712637.bam

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/038/SRR13712638/SRR13712638_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/038/SRR13712638/SRR13712638_2.fastq.gz
echo "Alineando..."
bowtie2 -x celegans -1 SRR13712638_1.fastq.gz -2 SRR13712638_2.fastq.gz -S SRR13712638.sam
rm SRR13712638_1.fastq.gz SRR13712638_2.fastq.gz
samtools view -bS SRR13712638.sam > SRR13712638.bam
rm SRR13712638.sam
samtools sort SRR13712638.bam -o SRR13712638.sorted.bam
rm SRR13712638.bam

rm celegans.fna celegans.b*

echo "Ejecutando script de R para construir la matriz de conteos"
Rscript script.r

if [ $? -eq 0 ]; then  # Si el script R terminó correctamente
    echo "Limpiando archivos BAM..."
    rm -f *sorted.bam *sorted.bam.bai
else
    echo "Error en el análisis R. Conservando archivos BAM para depuración."
fi


