
##FOR ONT
echo "Processing ONT data..."
for file in *_ONT.fastq
do
    ##ont use minimap for mapping
    minimap2 -t 40 -R "@RG\tID:"${file%_ONT.fastq}"\tSM:${file%_ONT.fastq}\tPL:nanopore" -a -x map-ont  tbdb.fasta $file | samtools sort -@ 20 -o ${file%_ONT.fastq}_ONT.bam
    samtools index -@ 50 ${file%_ONT.fastq}_ONT.bam
    ##variant calling
    freebayes -f tbdb.fasta --ploidy 1 --min-alternate-fraction 0.7 --haplotype-length -1 ${file%_ONT.fastq}_ONT.bam > ${file%_ONT.fastq}_ONT.vcf
    echo "First filters"
    ##filtering for low quality SNPs
    bcftools filter -i "FMT/DP>=10" ${file%_ONT.fastq}_ONT.vcf -Ov -o ${file%_ONT.fastq}_DPfilt.vcf
    bcftools +setGT ${file%_ONT.fastq}_DPfilt.vcf -- -i 'FMT/DP<5' -t q -n . > ${file%_ONT.fastq}_setGT.vcf
    bcftools view -i 'GT="1"' ${file%_ONT.fastq}_setGT.vcf -Ov -o ${file%_ONT.fastq}_GT1.vcf
    echo "Running bcftools to filter modlin"
    ##filter the regions
    bcftools view -T ^modlin.bed ${file%_ONT.fastq}_ONT.vcf > ${file%_ONT.fastq}_ONT_filtered.vcf
    cp ${file%_ONT.fastq}_ONT.vcf ${file%_ONT.fastq}_ONT_original.vcf
done

##FOR ILLUMINA
echo "Processing Illumina data..."
for file1 in *_ILL.1.fq
do
    file2=${file1%.1.fq}.2.fq

    
    trimmomatic PE -threads 20 -phred33 $file1 $file2 ${file1%.fq}_paired.fq ${file1%.fq}_unpaired.fq ${file2%.fq}_paired.fq ${file2%.fq}_unpaired.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    
     bwa mem -t 50 -K 10000000 -c 100 \
    -R "@RG\tID:${file1%_ILL.1.fq}\tSM:${file1%_ILL.1.fq}\tPL:illumina" -M -T 50 tbdb.fasta \
    ${file1%.fq}_paired.fq ${file2%.fq}_paired.fq | \
    samtools sort -n -@ 50 -o - | \
    samtools fixmate -@ 50 -m - - | \
    samtools sort -@ 50 -o - | \
    samtools markdup -@ 50 - ${file1%_ILL.1.fq}_ILL_marked.bam

    samtools index -@ 50 ${file1%_ILL.1.fq}_ILL_marked.bam

    echo "Running FreeBayes..."
    
    freebayes -f tbdb.fasta  --ploidy 1 --haplotype-length -1  ${file1%_ILL.1.fq}_ILL_marked.bam > ${file1%_ILL.1.fq}_ILL.vcf
    echo "First filters"
    bcftools filter -i "FMT/DP>=10" ${file1%_ILL.1.fq}_ILL.vcf -Ov -o ${file1%_ILL.1.fq}_DPfilt.vcf
    bcftools +setGT ${file1%_ILL.1.fq}_DPfilt.vcf -- -i 'FMT/DP<5' -t q -n . > ${file1%_ILL.1.fq}_setGT.vcf
    bcftools view -i 'GT="1"' ${file1%_ILL.1.fq}_setGT.vcf -Ov -o ${file1%_ILL.1.fq}_GT1.vcf
    echo "Running bcftools to filter modlin"
    bcftools view -T ^modlin.bed ${file1%_ILL.1.fq}_GT1.vcf > ${file1%_ILL.1.fq}_ILL_filtered.vcf
    cp ${file1%_ILL.1.fq}_GT1.vcf ${file1%_ILL.1.fq}_ILL_original.vcf
    
done

echo "Processing complete."
##############################ONT SVs#######################
#delly and sniffles structural variant calling
delly call -g tbdb.fasta S1_ONT.bam > S1_ONT.delly.vcf
sniffles -i S1_ONT.bam -v S1_ONT.sniffles.vcf

#filter delly and sniffles into beds to compare in IGV

bcftools view -i 'SVTYPE="DEL" && (END-POS)>15 && (END-POS)<50000 && GT="1/1"' S1_ONT.delly.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ONT.delly.bed
##for sniffles as well
bcftools view -i 'SVTYPE="DEL" && (END-POS)>15 && (END-POS)<50000 && GT="1/1"' S1_ONT.sniffles.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ONT.sniffles.bed

##############ILLUMINA SVs##############################

delly call -g tbdb.fasta S1_ILL_marked.bam > S1_ILL.delly.vcf
bcftools view -i 'SVTYPE="DEL" && (END-POS)>15 && (END-POS)<50000 && GT="1/1"' S1_ILL.delly.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ILL.delly.bed
