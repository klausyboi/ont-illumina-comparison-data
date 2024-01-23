
#############################################ONT###################################################
#mapping and calling
minimap2 -t 1 -R '@RG\tID:S1_ONT\tSM:S1_ONT\tPL:nanopore' -a -x map-ont tbdb.fasta S1_ONT.fastq | samtools sort -@ 50 -o S1_ONT.bam -
samtools index -@ 50 S1_ONT.bam
freebayes -f tbdb.fasta --haplotype-length -1 -F 0.7 S1_ONT.bam

#setgt to missing if the depth is lower than 5

bcftools view -T modlin.bed -s S1_ILL -x S1_ONT.vcf.gz | bcftools +setGT - -- -i 'FMT/DP<5' -t q -n .

#delly and sniffles structural variant calling
delly call -g tbdb.fasta S1_ONT.bam > S1_ONT.delly.vcf
sniffles -i S1_ONT.bam -v S1_ONT.sniffles.vcf
# filter delly and sniffles into beds to compare in IGV

bcftools view -c 2 S1_ONT.delly.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ONT.delly.bed

############################################ILLUMINA###############################################
bwa mem -t 1 -K 10000000 -C 100 -R '@RG\tID:S1_ILL\tSM:S1_ILL\tPL:illumina' -M -T 50 tbdb.fasta S1_ILL.1.fq S1_ILL.2.fq | samtools sort -@ 1 -o S1_ILL.bam -
samtools index -@ 50 S1_ILL.bam
freebayes -f tbdb.fasta --haplotype-length -1  S1_ILL.bam

bcftools view -T modlin.bed -s S1_ILL -x S1_ILL.vcf.gz | bcftools +setGT - -- -i 'FMT/DP<5' -t q -n .

delly call -g tbdb.fasta S1_ILL.bam > S1_ILL.delly.vcf

bcftools view -c 2 S1_ILL.delly.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ILL.delly.bed