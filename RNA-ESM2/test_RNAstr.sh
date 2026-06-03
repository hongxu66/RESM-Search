#######################################
# File Name: test_tmp.sh
# Author: hongxu
# mail: hongxu@szbl.ac.cn
#Created Time: Mon 15 Jan 2024 10:21:05 PM CST
###############################
#!/bin/bash

# filepath=/lustre/home/xhong/PostDotWorks/2024/MARS/database/DB/
filepath=/mnt/hdd/xhong/PostWork/2023/BPmap-LM/RNAstralign/RNAStrAlign/fasta/
ls $filepath | while read i
do
	ls $filepath$i | while read j;
	do
		# echo python esm_profile_test_batch.py produce.seqfile=${filepath}${i}/${j} produce.outname=${i}_${j%%.*}
		# python esm_profile_test_batch.py produce.seqfile=${filepath}${i}/${j} produce.outname=${j%%.*}
		echo python esm_profile_test_batch.py produce.seqfile=${filepath}${i}/${j} produce.outname=${j%%.*}
	done
		
done
