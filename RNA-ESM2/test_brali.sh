#######################################
# File Name: test_tmp.sh
# Author: hongxu
# mail: hongxu@szbl.ac.cn
#Created Time: Mon 15 Jan 2024 10:21:05 PM CST
###############################
#!/bin/bash

# filepath=/lustre/home/xhong/PostDotWorks/2024/MARS/database/DB/
filepath=/mnt/hdd/xhong/yikun/BRAliBase2/fasta/
ls $filepath | while read i
do
	python esm_profile_test_batch.py produce.seqfile=${filepath}${i} produce.outname=${i%%.*}
		
done
