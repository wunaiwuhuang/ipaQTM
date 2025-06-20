#! /bin/bash

cd /data1/wangwenhui/puQTM/1000_Genomes_phase3
for i in $(seq 1 22)
do
	echo $i
	
	pt="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chr${i}.GRCh38.GT.crossmap.vcf.gz"
	wget ${pt}

done

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chrX.GRCh38.GT.crossmap.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chrY.GRCh38.GT.crossmap.vcf.gz



