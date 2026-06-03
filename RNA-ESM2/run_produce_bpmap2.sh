source ../LModel/bin/activate
# for((i=0;i<81;i++));do python esm2embatt_650M_con.py produce.seqfile=/mnt/hdd/xhong/yikun/RNAcentral/dir${i}/seq.fa produce.outname=dir${i}_out;done
# python esm2embatt_650M_con.py produce.seqfile=/mnt/hdd/xhong/yikun/RNAcentral/dir23/seq.fa produce.outname=dir23_out
# CUDA_VISIBLE_DEVICES=3 python esm2embatt_650M_con.py produce.seqfile=/mnt/hdd/xhong/yikun/RNAcentral/dir74/seq.fa produce.outname=dir74_out
# CUDA_VISIBLE_DEVICES=2 python esm2embatt_650M_con.py produce.seqfile=/mnt/hdd/xhong/yikun/RNAcentral/dir62/seq.fa produce.outname=dir62_out
# CUDA_VISIBLE_DEVICES=1 python esm2embatt_650M_con.py produce.seqfile=/mnt/hdd/xhong/yikun/RNA-ESM2/all_seq.fa produce.outname=benchmark_out

input=$1

dirpath=/mnt/remote_home/xhong/PostWork/2025/latest_pos/$input

for seqfile in ${dirpath}/*.fa
do
    id=${seqfile##*/}
    if [ ! -f /mnt/remote_home/xhong/PostWork/2024/RNA-ESM2/results_650_con/${input}_${id}_allmaps_reverseseq.bpmap ];
    then
    # echo CUDA_VISIBLE_DEVICES=1 python esm2embatt_650M_con.py produce.seqfile=${seqfile} produce.outname=chr1_${id}
        CUDA_VISIBLE_DEVICES=4 python esm2embatt_650M_con_reverseseq.py produce.seqfile=${seqfile} produce.outname=${input}_${id}
        # CUDA_VISIBLE_DEVICES=3 python esm2embatt_650M_con.py produce.seqfile=${seqfile} produce.outname=${input}_${id}
        # CUDA_VISIBLE_DEVICES=3 python esm2embatt_650M_con.py produce.seqfile=${seqfile} produce.outname=${input}_${id}
    fi
done

