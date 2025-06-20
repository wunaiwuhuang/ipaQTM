#! /bin/bash 


cd /data1/wangwenhui/puQTM/1000_Genomes_phase3
for i in $(seq 1 22)
do
        echo $i

	vcf_file="/data1/wangwenhui/puQTM/1000_Genomes_phase3/phase3.chr${i}.GRCh38.GT.crossmap.vcf.gz"
	plink --vcf ${vcf_file} --freq --memory 50000 --out /data1/wangwenhui/puQTM/1000_Genomes_phase3/1000GP.hg38.chr${i}

done

