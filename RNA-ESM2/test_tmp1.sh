#######################################
# File Name: test_tmp.sh
# Author: hongxu
# mail: hongxu@szbl.ac.cn
#Created Time: Mon 15 Jan 2024 10:21:05 PM CST
###############################
#!/bin/bash

# filepath=/lustre/home/xhong/PostDotWorks/2024/MARS/database/DB/
filepath=/home/xhong/MARS/DB/
echo `printf "%03d\n" $1`
printf "%03d\n" $1 | while read i;
do
	ls -l ${filepath}rdp2mdb5_${i}.fasta.split/* | awk '$5>200000000{print $NF}' | while read j;
	do 
		python esm_profile_test_batch.py produce.seqfile=${j} produce.outname=${j##*/}
		
	done
done
