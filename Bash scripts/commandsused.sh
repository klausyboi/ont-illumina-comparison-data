
#mapping and calling
 minimap2 -t 1 -R '@RG\tID:tbprofiler\tSM:tbprofiler\tPL:nanopore' -a -x map-ont  S1_ONT.fastq | samtools sort -@ 1 -o S1_ONT.bam
 samtools index -@ 1 S1_ONT.bam
 freebayes -f tbdb.fasta --haplotype-length -1 -F 0.7 S1_ONT.bam

#setgt to missing if the depth is lower than 5

bcftools view -T modlin.bed  S1_ONT.vcf.gz -s S1_ONT,S1_ILL -x | bcftools +setGT - -- -i 'FMT/DP<5' -t q -n . 

#delly and sniffles structural variant calling
delly call -g tbdb.fasta S1_ONT.bam > S1_ONT.delly.vcf
sniffles -i S1_ONT.bam -v S1_ONT.sniffles.vcf
# filter delly and sniffles into beds to compare in IGV

bcftools view -c 2 S1_ONT.delly.vcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > S1_ONT.delly.bed

