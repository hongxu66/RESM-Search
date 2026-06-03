import os,sys
import re
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
    model_path : str = str(current_directory / "ckpt/RNA-ESM2-trans-2a100-mappro-KDNY-epoch=07-valid_F1=0.561-v1.ckpt")
    device: str = "cpu"
    # device: str = "gpu"

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
    embed_dim: int = 640
    num_attention_heads: int = 20
    dropout: float = 0.1
    attention_dropout: float = 0.1
    activation_dropout: float = 0.1
    attention_type: str = "standard"
    performer_attention_features: int = 256
    num_layers: int = 30
    max_seqlen: int = 1024

@dataclass
class TrainConfig:
    pass

@dataclass
class LoggingConfig:
    pass

@dataclass
class Config:
    data: DataConfig = DataConfig()
    train: TrainConfig = TrainConfig()
    model: TransformerConfig = TransformerConfig()
    optimizer: OptimizerConfig = OptimizerConfig()
    logging: LoggingConfig = LoggingConfig()

cs = ConfigStore.instance()
cs.store(name="config", node=Config)
cs.store(group="data", name="default", node=DataConfig)
cs.store(group="train", name="default", node=TrainConfig)
cs.store(group="optimizer", name="default", node=OptimizerConfig)
cs.store(group="model", name="emb-transformer", node=TransformerConfig)
cs.store(group="logging", name="default", node=LoggingConfig)

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
    def readfasta1(filename):
        seqs = []
        tem = bseq.parse(filename,'fasta')
        for item in tem:
            seqs.append((str(item.id),str(item.seq)))
        return seqs

    model_path = cfg.data.model_path
    device = torch.device(cfg.data.device)
    _, protein_alphabet = rna_esm.pretrained.esm2_t30_150M_UR50D()
    rna_alphabet = rna_esm.data.Alphabet.from_architecture(cfg.data.architecture)

    protein_vocab = Vocab.from_esm_alphabet(protein_alphabet)
    rna_vocab = Vocab.from_esm_alphabet(rna_alphabet)
    rna_map_dict = mapdict(protein_vocab, rna_vocab)
    rna_map_vocab = Vocab.from_esm_alphabet(rna_alphabet, rna_map_dict)

    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/yikun/RNA-ESM2/testset_allseqs.fasta')
    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/PDBSS/fasta/output.fasta')
    # seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/compare2RNAcmap3/benchmark/all_sequences.fasta')
    seqs = readfasta1('/lustre/home/xhong/PostDotWorks/2023/compare2RNAcmap3/src/dir1/seq.fa')
    # seqs = readfasta1(sys.argv[1])
#    seq = "GGGGACUCCAGAGGUCGAGAGACCGGAGAUAUCACCC"
#    tokens = torch.from_numpy(vocab.encode(seq))

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

    datas = []
    conts = []
    attens = []
    with torch.no_grad():
        for item in seqs:
            id,seq = item
            print(id,seq)
            tokens = torch.from_numpy(rna_map_vocab.encode(seq))
            tokens = tokens.unsqueeze(0)
            logits = model(tokens.to(device))["logits"]
            start_idx = int(rna_map_vocab.prepend_bos) ##
            end_idx = logits.size(-2) - int(rna_map_vocab.append_eos) ##
#            logits = logits[:,start_idx:end_idx,:]
            print(logits.shape)
            outfea = logits[:,start_idx:end_idx,4:8] ##
#            outfea = logits[:,:,4:8] ##
            outfea = outfea.squeeze()
#            print(outfea)
            outfea = F.softmax(outfea,dim=1)
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
            print(predictions.shape)
            conts.append((id,seq,predictions))
            # attens.append((id,seq,predictions1))
    np.savez('./benchmark_features.npz',a=datas)
    # np.savez('./benchmark_attentions.npz',a=attens)
    np.savez('./benchmark_contacts.npz',a=conts)

@hydra.main(config_name="config")
def produce_seqprf_rfam(cfg: Config) -> None:
    def readfasta1(filename):
        seqs = []
        tem = bseq.parse(filename,'fasta')
        for item in tem:
            # for base in str(item.seq):
            # if len(list(set(['A','U','C','G','N','T']) | set(list(str(item.seq))))) > 6:
            #     continue
            seqs.append((str(item.id),re.sub(r'[^AUCGTN]','N',str(item.seq))))
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
    # produce_seqprf()
    produce_seqprf_rfam()


