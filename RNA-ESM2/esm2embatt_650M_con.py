import os,sys
import re
import math
import json
from typing import Tuple
from pathlib import Path
from pytorch_lightning import seed_everything
import torch
import rna_esm
from evo.tokenization import Vocab, mapdict
from model import ESM2
from dataclasses import dataclass
import hydra
from hydra.core.config_store import ConfigStore
from omegaconf import DictConfig, OmegaConf

import Bio.SeqIO as bseq
import numpy as np
import torch.nn.functional as F
from sklearn.preprocessing import normalize


seed_everything(42)

current_directory = Path(__file__).parent.absolute()


@dataclass
class DataConfig:
    architecture: str = "rna-esm"
    num_workers: int = 16
    # model_path : str = str(current_directory / "ckpt/RNA-ESM2-trans-2a100-mappro-KDNY-epoch=07-valid_F1=0.561-v1.ckpt")
    # model_path : str = str(current_directory / "ckpt/RNA-ESM2-trans-2a100-mappro-KDNY-epoch=06-valid_F1=0.564.ckpt")
    model_path : str = str(current_directory / "ckpt/RUNGPURNA-ESM2-trans-mappro-KDNY-epoch-30-step-15650.pt")
    # device: str = "cpu"
    device: str = "cuda:0"

@dataclass
class OptimizerConfig:
    name: str = "adam"
    learning_rate: float = 4e-4
    weight_decay: float = 1e-2     #3e-4
    lr_scheduler: str = "warmup_linear" #"warmup_cosine"
    warmup_steps: int = 16000    #16000
    adam_betas: Tuple[float, float] = (0.9, 0.98)  # 0.9, 0.999
    max_steps: int = 10000000

@dataclass
class TransformerConfig:
    embed_dim: int = 1280
    num_attention_heads: int = 20
    dropout: float = 0.1
    attention_dropout: float = 0.1
    activation_dropout: float = 0.1
    attention_type: str = "standard"
    performer_attention_features: int = 256
    num_layers: int = 33
    max_seqlen: int = 1024

@dataclass
class TrainConfig:
    pass

@dataclass
class LoggingConfig:
    pass

@dataclass
class ProduceConfig:
    seqfile: str = "All_IRES.fa"
    outname: str = "All_IRES"

@dataclass
class Config:
    data: DataConfig = DataConfig()
    train: TrainConfig = TrainConfig()
    model: TransformerConfig = TransformerConfig()
    optimizer: OptimizerConfig = OptimizerConfig()
    logging: LoggingConfig = LoggingConfig()
    produce: ProduceConfig = ProduceConfig()

cs = ConfigStore.instance()
cs.store(name="config", node=Config)
cs.store(group="data", name="default", node=DataConfig)
cs.store(group="train", name="default", node=TrainConfig)
cs.store(group="optimizer", name="default", node=OptimizerConfig)
cs.store(group="model", name="emb-transformer", node=TransformerConfig)
cs.store(group="logging", name="default", node=LoggingConfig)
cs.store(group="produce",name="produce",node=ProduceConfig)

@hydra.main(config_name="config")
def TS_test(cfg: Config) -> None:
    model_path = cfg.data.model_path
    device = torch.device(cfg.data.device)
    _, protein_alphabet = rna_esm.pretrained.esm2_t30_150M_UR50D()
    rna_alphabet = rna_esm.data.Alphabet.from_architecture(cfg.data.architecture)

    protein_vocab = Vocab.from_esm_alphabet(protein_alphabet)
    rna_vocab = Vocab.from_esm_alphabet(rna_alphabet)
    rna_map_dict = mapdict(protein_vocab, rna_vocab)
    rna_map_vocab = Vocab.from_esm_alphabet(rna_alphabet, rna_map_dict)

    # seq = "GGGGACUCCAGAGGUCGAGAGACCGGAGAUAUCACCC"
    seq = ["GGGGACUCCAGAGGUCGAGAGACCGGAGAUAUCACCC","GGGGAAAAAAUUUUUCCAA"]
    tokens = torch.from_numpy(rna_map_vocab.encode(seq))

    model = ESM2(
        vocab=protein_vocab,
        model_config=cfg.model,
        optimizer_config=cfg.optimizer,
        contact_train_data=None,
        token_dropout=True,
    )

    model.load_state_dict(torch.load(
        model_path,
        map_location=device)['state_dict'], strict=True)
    model = model.eval()
    model = model.to(device)

    with torch.no_grad():
        tokens = tokens.unsqueeze(0)
        # logits = model(tokens.to(device))["logits"]
        logits = model(tokens.to(device))
        print(logits.keys())
        aaaa = logits['representations']
        print('aaaa',aaaa)
        logits = logits["logits"]
        start_idx = int(rna_map_vocab.prepend_bos)
        end_idx = logits.size(-2) - int(rna_map_vocab.append_eos)
        logits = logits[:, start_idx:end_idx, :]
        predictions = model(tokens.to(device),return_contacts=True)['attentions']
        predictions = predictions[:,:,:,start_idx:end_idx,start_idx:end_idx]
        predictions = predictions.squeeze(0)
        predictions = predictions.numpy()
        predictions1 = model(tokens.to(device),return_contacts=True)["contacts"]
        # predictions1 = predictions1[:,start_idx:end_idx,start_idx:end_idx]
        predictions1 = predictions1[:,:,:]
		
        print(logits.shape)
        print(logits)
        print(len(seq))
        print(predictions.shape)
        print(predictions1.shape)

@hydra.main(config_name="config")
def produce_seqprf(cfg: Config) -> None:
    def load_model_without_module(model_state_file):
        from collections import OrderedDict
        new_checkpoint = OrderedDict()
        checkpoint = torch.load(model_state_file, map_location="cpu")

        for k, v in checkpoint.items():
            name = k[7:] # remove module.
            new_checkpoint[name] = v
        return new_checkpoint

    def writejson(data, outfile):
        with open(outfile, 'w') as f:
            json.dump(data, f)

    def readjson(jsonfile):
        with open(jsonfile) as f:
            data = json.load(f)
        return data

    def readfasta2(filename):
        seqs = []
        tem = bseq.parse(filename,'fasta')
        outname = 'seq'
        for i,item in enumerate(tem):
            seqs.append((outname + str(i),re.sub(r'[^AUCGTN]','N',str(item.seq))))
        return seqs
    

    def readfasta1(filename):
        seqs = []
        tem = bseq.parse(filename,'fasta')
        for item in tem:
            if len(item.seq) > 1024:
                continue
            # seqs.append((str(item.id),re.sub(r'[^AUCGTN]','N',str(item.seq).upper()[:70])))
            seqs.append((str(item.id),re.sub(r'[^AUCGTN]','N',str(item.seq).upper())))
        return seqs

    def splitseq(seqs):
        ''' the max sequence length of language model required is 1024'''
        maxLength = 1024
        cutLength = 1000
        overlap = 200
        mapdict = {}
        seqdict = {}
        newseqs = []
        index = 0
        for item in seqs:
            id = 'seq_' + str(index)
            if len(item[1]) > maxLength:
                num = math.ceil((len(item[1]) - overlap)/ (cutLength - overlap))
                for j in range(num):
                    id1 = id + '-' + str(j)
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

    def pair2map(pairs):
        s = ''
        for item in pairs:
            s += ' '.join(['CON'] + [str(i) for i in item])
            s += '\n'
        return s

    def prf2map(seq, prf, ss=None):
        prf = prf.T
        base = 'AUCG'
        index1 = [0, 3, 2, 1]
        s = ''
        index = 0
        for i in range(len(seq)):
            # pA, pU, pC, pG = (1 - init) / 3, (1 - init) / 3, (1 - init) / 3, (1 - init) / 3
            if seq[i] == '-':
                continue
            if ss == None:
                s += ' '.join(['PRF', str(index), seq[i], 'N'] + [str(prf[j, i]) for j in index1] + ['\n'])
            else:
                if ss[i] == '.':
                    s += ' '.join(['PRF', str(index), seq[i], 'N'] + [str(prf[j, i]) for j in index1] + ['\n'])
                else:
                    s += ' '.join(['PRF', str(index), seq[i], 'P'] + [str(prf[j, i]) for j in index1] + ['\n'])
            index += 1
        return s

    def outtargetMapLast(seq, prf, con,outname, cutoff=0.01):
        con = np.squeeze(con)
        newmatrix = con
        pairs = []
        L = len(seq)
        for i in range(L):
            for j in range(i + 2, L):
                if newmatrix[i, j] > cutoff:
                    pairs.append((i, j, newmatrix[i, j]))
        s = '##' + outname + '\n'
        # s = ''
        s += 'LEN {}\n'.format(len(seq))
        s += pair2map(pairs)
        s += prf2map(seq, prf)
        # tem = open(outname, 'w+')
        # tem.write(s)
        # tem.close()
        return s


    print(OmegaConf.to_yaml(cfg))
    print('cuda:::',torch.cuda.is_available())
    outname = cfg.produce.outname

    model_path = cfg.data.model_path
    device = torch.device(cfg.data.device)
    # _, protein_alphabet = rna_esm.pretrained.esm2_t30_150M_UR50D()
    _, protein_alphabet = rna_esm.pretrained.esm2_t33_650M_UR50D()
    rna_alphabet = rna_esm.data.Alphabet.from_architecture(cfg.data.architecture)
    
    # kk = ProduceConfig(sys.argv[1],sys.argv[2])

    protein_vocab = Vocab.from_esm_alphabet(protein_alphabet)
    rna_vocab = Vocab.from_esm_alphabet(rna_alphabet)
    rna_map_dict = mapdict(protein_vocab, rna_vocab)
    rna_map_vocab = Vocab.from_esm_alphabet(rna_alphabet, rna_map_dict)

    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/yikun/RNA-ESM2/testset_allseqs.fasta')
    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/PDBSS/fasta/output.fasta')
    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/compare2RNAcmap3/benchmark/all_sequences.fasta')
    # allseqs = readfasta1(cfg.produce.seqfile)
    # print(allseqs)
    # seqmapdict,seqdict,seqs = splitseq(allseqs)
    # print(seqs)
    # writejson(seqmapdict,outname + '_mapdict.json')
    # writejson(seqdict,outname + '_seqdict.json')
    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/compare2RNAcmap3/src/dir1/seq.fa')
    # seqs = readfasta1(sys.argv[1])
    allseqs = readfasta1(cfg.produce.seqfile)
    # allseqs = readfasta2(cfg.produce.seqfile)
    #### 
    # seqmapdict,seqdict,seqs = splitseq(allseqs)
    
#    seq = "GGGGACUCCAGAGGUCGAGAGACCGGAGAUAUCACCC"
#    tokens = torch.from_numpy(vocab.encode(seq))

    model = ESM2(
        vocab=protein_vocab,
        model_config=cfg.model,
        optimizer_config=cfg.optimizer,
        contact_train_data=None,
        token_dropout=True,
    )


    model.load_state_dict(load_model_without_module(model_path))

    model = model.eval()
    model = model.to(device)

    datas = []
    conts = []
    attens = []
    allmaps = ''
    # results_path = '/mnt/hdd/xhong/yikun/RNA-ESM2/results/'
    # results_path = '/mnt/remote_home/xhong/PostWork/2024/RNA-ESM2/results_650_con/'
    # results_path = '/mnt/remote_home/xhong/PostWork/2026/human_genome_contact_map/'
    results_path = '/mnt/remote_home/xhong/PostWork/2025/flav_genomes/map_80/'
    with torch.no_grad():
        # for item in seqs:
        for item in allseqs:
            id,seq = item
            print(id,seq)
            if len(seq) <= 1:
                continue
            # embedding_path = os.path.join(results_path,id + '_emb.npy')
            # attentions_path = os.path.join(results_path,id + '_atp.npy')
            # if os.path.exists(embedding_path):
            #     continue
            tokens = torch.from_numpy(rna_map_vocab.encode(seq))
            tokens = tokens.unsqueeze(0)
            # logits = model(tokens.to(device))["logits"]
            outs = model(tokens.to(device),repr_layers=[33],return_contacts=True)
            logits = outs['logits']
            start_idx = int(rna_map_vocab.prepend_bos) ##
            end_idx = logits.size(-2) - int(rna_map_vocab.append_eos) ##
            # emb = outs['representations'][33][:,start_idx:end_idx,:]
            # emb = emb.squeeze()
            # try:
            #     emb = F.softmax(emb,dim=1)
            # except IndexError:
            #     pass
            # emb = emb.cpu().numpy()
            # att = outs['attentions'].squeeze(0)[:,:,start_idx:end_idx,start_idx:end_idx].permute(2,3,0,1).cpu().numpy()
            # np.save(embedding_path,emb)
            # np.save(attentions_path,att)
#             print(logits.shape)
            outfea = logits[:,start_idx:end_idx,4:8] ##
#            outfea = logits[:,:,4:8] ##
            outfea = outfea.squeeze()
#            print(outfea)
            outfea = F.softmax(outfea,dim=1)
            outfea = outfea.cpu()
            outfea = outfea.numpy()
            # datas.append((id,seq,outfea))
            # predictions = model(tokens.to(device),return_contacts=True)["contacts"]
            predictions = outs['contacts']
            predictions = predictions[:,:,:]
            predictions = predictions.squeeze(0) 
            predictions = predictions.cpu() 
            predictions = predictions.numpy()
            allmaps += outtargetMapLast(seq,outfea,predictions,id)

            conts.append((id,seq,predictions))

    # np.savez('./benchmark_features.npz',a=datas)
    # np.savez('./' + outname + '_features.npz',a=datas)
    # np.savez('./benchmark_attentions.npz',a=attens)
    # np.savez('./benchmark_contacts.npz',a=conts)
    # np.savez('./' + outname + '_contacts.npz',a=conts)
    # return 0
    temout = open(results_path + outname + '_allmaps.bpmap','w+')
    temout.write(allmaps)
    temout.close()
    # return 0
    # writejson(seqmapdict,results_path + '/' + outname + '_mapdict.json')
    # writejson(seqdict,results_path + '/' + outname + '_seqdict.json')

def load_model_hwei(rna_esm2_model_path, prot_esm2_model_path):

    def load_model_without_module(model_state_file):
        from collections import OrderedDict
        new_checkpoint = OrderedDict()
        checkpoint = torch.load(model_state_file, map_location="cpu")

        for k, v in checkpoint.items():
            name = k[7:] # remove module.
            new_checkpoint[name] = v
        return new_checkpoint
    

    model, _ = esm.pretrained.load_model_and_alphabet_local(prot_esm2_model_path)
    model.load_state_dict(load_model_without_module(rna_esm2_model_path))

    return model


@hydra.main(config_name="config")
def produce_seqprf_rfam(cfg: Config) -> None:
    def readfasta1(filename):
        seqs = []
        tem = bseq.parse(filename,'fasta')
        for item in tem:
            # for base in str(item.seq):
            # if len(list(set(['A','U','C','G','N','T']) | set(list(str(item.seq))))) > 6:
            #     continue
            seqs.append((str(item.id),re.sub(r'[^AUCGTN]','N',str(item.seq).upper())))
			# re.sub(r'[^AUCGTN]','N',a)
        return seqs

    def readseqdir():
        dirpath='/lustre/home/xhong/PostDotWorks/2023/compare2RNAcmap3/src/'
        seqdict = {}
        for item in os.listdir(dirpath):
            if 'dir' in item:
                seqdict[item] = readfasta1(dirpath + item + '/seq.fa')
        return seqdict

    def readpathother(dirpath):
        seqdict = {}
        dirlists = os.listdir(dirpath)
        for item in dirlists:
            if os.path.exists('/lustre/home/xhong/PostDotWorks/2023/BPmap_new_training/LM/' + item + '_contacts.npz'):
                continue
            if 'seq' in item:
                seqdict[item] = readfasta1(dirpath + item)
        return seqdict

    def readpath(dirpath):
        seqdict = {}
        dirlists = os.listdir(dirpath)
        for item in dirlists:
            if os.path.exists('/lustre/home/xhong/PostDotWorks/2023/BPmap_new_training/LM/' + item + '_contacts.npz'):
                continue
            seqdict[item] = readfasta1(dirpath + item + '/allseqs.fasta')
        return seqdict

    model_path = cfg.data.model_path
    device = torch.device(cfg.data.device)
    _, protein_alphabet = rna_esm.pretrained.esm2_t30_150M_UR50D()
    rna_alphabet = rna_esm.data.Alphabet.from_architecture(cfg.data.architecture)

    protein_vocab = Vocab.from_esm_alphabet(protein_alphabet)
    rna_vocab = Vocab.from_esm_alphabet(rna_alphabet)
    rna_map_dict = mapdict(protein_vocab, rna_vocab)
    rna_map_vocab = Vocab.from_esm_alphabet(rna_alphabet, rna_map_dict)

    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/yikun/RNA-ESM2/testset_allseqs.fasta')
    # seqs = readfasta1(sys.argv[1])
#    seq = "GGGGACUCCAGAGGUCGAGAGACCGGAGAUAUCACCC"
#    tokens = torch.from_numpy(vocab.encode(seq))
    # seqdict = readpath('/lustre/home/xhong/PostDotWorks/2022/RFAM/fasta/')
    # seqdict = readpathother('/lustre/home/xhong/PostDotWorks/2023/PDBSS/fasta/')
    seqdict = readseqdir()

    model = ESM2(
        vocab=protein_vocab,
        model_config=cfg.model,
        optimizer_config=cfg.optimizer,
        contact_train_data=None,
        token_dropout=True,
    )

    model.load_state_dict(torch.load(
        model_path,
        map_location=device)['state_dict'], strict=True)
    model = model.eval()
    model = model.to(device)


    with torch.no_grad():
        for keyword in seqdict.keys():
            datas = []
            conts = []
            attens = []
            print(keyword)
            for item in seqdict[keyword]:
                id,seq = item
                tokens = torch.from_numpy(rna_map_vocab.encode(seq))
                tokens = tokens.unsqueeze(0)
                logits = model(tokens.to(device))["logits"]
                start_idx = int(rna_map_vocab.prepend_bos) ##
                end_idx = logits.size(-2) - int(rna_map_vocab.append_eos) ##
    #            logits = logits[:,start_idx:end_idx,:]
                # print(logits.shape)
                outfea = logits[:,start_idx:end_idx,4:8] ##
    #            outfea = logits[:,:,4:8] ##
                outfea = outfea.squeeze()
    #            print(outfea)
                try:
                    outfea = F.softmax(outfea,dim=1)
                except IndexError:
                    continue
    #            outfea = normalize(outfea,axis=1,norm='l1')
                outfea = outfea.numpy()
                # print(outfea)
                datas.append((id,seq,outfea))
                predictions = model(tokens.to(device),return_contacts=True)["contacts"]
                # predictions = predictions[:,start_idx:end_idx,start_idx:end_idx]
                predictions = predictions[:,:,:]
                predictions = predictions.squeeze(0)
                predictions = predictions.numpy()
                # predictions1 = model(tokens.to(device),return_contacts=True)["attentions"]
                # predictions1 = predictions1[:,:,:,start_idx:end_idx,start_idx:end_idx]
                # predictions1 = predictions1.squeeze(0)
                # predictions1 = predictions1.numpy()
                # print(predictions.shape)
                conts.append((id,seq,predictions))
                # attens.append((id,seq,predictions1))
            np.savez('./' + keyword + '_features.npz',a=datas)
            # np.savez('./' + keyword + '_attentions.npz',a=attens)
            np.savez('./' + keyword + '_contacts.npz',a=conts)


if __name__ == "__main__":
    # TS_test()
    produce_seqprf()
    # produce_seqprf_rfam()


