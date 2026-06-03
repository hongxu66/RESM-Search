
input=$1_neg

ls /mnt/remote_home/xhong/PostWork/2025/latest_pos/chr${input}_60/*.fa | while read i;
do 
    echo ${i##*/};
    CUDA_VISIBLE_DEVICES=4 python esm2embatt_650M_con.py produce.seqfile=/mnt/remote_home/xhong/PostWork/2025/latest_pos/chr${input}_60/${i##*/} produce.outname=chr${input}_60/${i##*/};
done
