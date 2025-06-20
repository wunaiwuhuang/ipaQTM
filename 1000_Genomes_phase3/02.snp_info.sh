#! /bin/bash


cd /data1/wangwenhui/puQTM/1000_Genomes_phase3
for i in $(seq 1 22)
do
        echo $i

        vcf_file="/data1/wangwenhui/puQTM/1000_Genomes_phase3/phase3.chr${i}.GRCh38.GT.crossmap.vcf.gz"
	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${vcf_file} > SNP.info.chr${i}.txt

done
