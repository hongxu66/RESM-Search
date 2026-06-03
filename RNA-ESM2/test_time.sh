date
# python esm_profile_test_batch.py produce.seqfile=/home/xhong/MARS/rdp2mdb5_149.fasta.split/rdp2mdb5_149.part_005.fasta produce.outname=rdp2mdb5_149.part_005
# python esm_profile_test_batch.py produce.seqfile=/home/xhong/MARS/DB/rdp2mdb5_001.fasta.split/rdp2mdb5_001.part_024.fasta produce.outname=rdp2mdb5_001.part_024
for i in `seq 5 48`;do printf "%03d\n" $i;done | while read j;
do 
    echo $j
    # python esm_profile_test_batch.py produce.seqfile=/home/xhong/MARS/DB/rdp2mdb5_147.fasta.split/rdp2mdb5_147.part_$j.fasta produce.outname=rdp2mdb5_147.part_$j
    # python esm_profile_test_batch.py produce.seqfile=/home/xhong/MARS/DB/rdp2mdb5_146.fasta.split/rdp2mdb5_146.part_$j.fasta produce.outname=rdp2mdb5_146.part_$j
    python esm_profile_test_batch.py produce.seqfile=/home/xhong/MARS/DB/rdp2mdb5_148.fasta.split/rdp2mdb5_148.part_$j.fasta produce.outname=rdp2mdb5_148.part_$j
done 
date
