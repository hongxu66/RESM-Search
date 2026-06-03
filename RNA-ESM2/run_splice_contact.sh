
input=$1

ls /mnt/remote_home/xhong/PostWork/2025/latest_pos/chr${input}_60/*.fa | while read i;
do 
    if [ -f /mnt/remote_home/xhong/PostWork/2026/human_genome_contact_map/chr${input}_60/${i##*/}_allmaps.bpmap ];
    then
        continue
    fi
    echo ${i##*/};
    CUDA_VISIBLE_DEVICES=2 python esm2embatt_650M_con.py produce.seqfile=/mnt/remote_home/xhong/PostWork/2025/latest_pos/chr${input}_60/${i##*/} produce.outname=chr${input}_60/${i##*/};
done
