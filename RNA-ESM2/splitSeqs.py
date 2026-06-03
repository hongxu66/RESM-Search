import os,sys
import math
import re
import numpy as np
import Bio.SeqIO as bseq
import json

def writejson(data, outfile):
    with open(outfile, 'w') as f:
        json.dump(data, f)

def readjson(jsonfile):
    with open(jsonfile) as f:
        data = json.load(f)
    return data


def readfasta1(filename):
    seqs = []
    tem = bseq.parse(filename,'fasta')
    for item in tem:
        # for base in str(item.seq):
        # if len(list(set(['A','U','C','G','N','T']) | set(list(str(item.seq))))) > 6:
        #     continue
        seqs.append((str(item.id),re.sub(r'[^AUCGTN]','N',str(item.seq).upper())))
        # seqs.append((str(item.id),str(item.seq)))
    		# re.sub(r'[^AUCGTN]','N',a)
    return seqs

def reverse_seq(seq):
    seq = list(seq)
    seq.reverse()
    tmpseq = ''
    for base in seq:
        if base == 'A':
            tmpseq += 'T'
        elif base == 'T' or base == 'U':
            tmpseq += 'A'
        elif base == 'C':
            tmpseq += 'G'
        elif base == 'G':
            tmpseq += 'C'
        else:
            tmpseq += 'N'
    return tmpseq

def splitseq(seqs):
    ''' the max sequence length of language model required is 1024'''
    maxLength = 128
    cutLength = 128
    overlap = 90
    mapdict = {}
    seqdict = {}
    newseqs = []
    index = 0
    for item in seqs:
        id = 'seq_' + str(index)
        if len(item[1]) > maxLength:
            num = math.ceil((len(item[1]) - overlap)/ (cutLength - overlap))
            for j in range(num):
                id1 = id + '_' + str(j)
                seq = item[1][j * (cutLength-overlap):(j + 1) * (cutLength - overlap) + overlap]
                mapdict[id1] = item[0]
                seqdict[id1] = seq
                newseqs.append((id1,seq))
            # id1 = id + '-' + str(num)
            # seq = item[1][num * (cutLength-overlap):]
            # mapdict[id1] = item[0]
            # seqdict[id1] = seq
            # newseqs.append((id1,seq))
        else:
            id += '-0'
            mapdict[id] = item[0]
            seqdict[id] = item[1]
            newseqs.append((id, item[1]))
        index += 1
    return mapdict, seqdict, newseqs

def main():
    seqfile = sys.argv[1]

    id = seqfile.split('.')[0]
    if not os.path.exists(id):
        os.mkdir(id)
    seq = readfasta1(seqfile)

    mapdict,seqdict,newseqs = splitseq(seq)
    writejson(mapdict,os.path.join(id,'mapdict.json'))
    writejson(seqdict,os.path.join(id,'seqdict.json'))

    index = -1
    tmpout = open(id + '/tmp' + str(index) + '.fa','w+')
    for i,item in enumerate(newseqs):
        if i % 2000 == 0:
            tmpout.close()
            index += 1
            tmpout = open(id + '/tmp' + str(index) + '.fa','w+')
        tmpid,tmpseq = item
        if tmpseq.count('N') == len(tmpseq):
            continue
        tmpout.write(f'>{tmpid}\n{tmpseq}\n')
    tmpout.close()

def test():
    seqfile = sys.argv[1]
    id = seqfile.split('.')[0]
    seq = readfasta1(seqfile)
    startindex = int(sys.argv[2])-1
    endindex = int(sys.argv[3])
    tmpseq = seq[0][1][startindex:endindex]
    if sys.argv[4] == '-':
        tmpseq = reverse_seq(tmpseq)
    print(f'>{id}:{sys.argv[2]}-{sys.argv[3]}\n{tmpseq}')
    
def test1():
    seqfile = sys.argv[1]
    seq = readfasta1(seqfile)
    startindex = int(int(int(sys.argv[2]) / 38) * 38) 
    endindex = startindex + 128
    tmpseq = seq[0][1][startindex:endindex]
    # print(startindex,endindex)
    # print(f'>test\n{tmpseq}')
    print(f'>test:{sys.argv[2]}\n{tmpseq}\n')

def test2():
    seqfile = sys.argv[1]
    seq = readfasta1(seqfile)
    indexlists = [int(i.rstrip())-1 for i in open(sys.argv[2]).readlines()]
    seqdict = {}
    for index in indexlists:
        code = int(index/38)
        startindex = int(int(index/38) * 38)
        endindex = startindex + 128
        code_seq = seq[0][1][startindex:endindex]
        # print(startindex,endindex,code_seq)
        if code not in seqdict.keys():
            seqdict[code] = [code_seq,1]
        else:
            num = seqdict[code][1] + 1 
            seqdict[code] = [code_seq,num]
    for code in seqdict.keys():
        code_seq,num = seqdict[code]
        print(f'>test:{code}:{num}\n{code_seq}')

def test3():
    seqfile = sys.argv[1]
    seq = readfasta1(seqfile)
    nlen = math.ceil((len(seq[0][1])-90)/38)
    print(nlen,len(seq[0][1]))
    # indexlists = [int(i.rstrip())-1 for i in open(sys.argv[2]).readlines()]
    seqdict = {}
    for index in range(nlen):
        code = index * 38
        # startindex = int(int(index/38) * 38)
        startindex = index * 38
        endindex = startindex + 128
        if endindex > len(seq[0][1]):
            endindex = len(seq[0][1])
        # print(startindex,endindex)
        code_seq = seq[0][1][startindex:endindex]
        # print(startindex,endindex,code_seq)
        if code not in seqdict.keys():
            seqdict[code] = code_seq
    for code in seqdict.keys():
        code_seq = seqdict[code]
        print(f'>test:{code}\n{code_seq}')

if __name__ == '__main__':
    # main()
    # test()
    # test1()
    # test2()
    test3()
