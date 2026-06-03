# RESM-Search

RESM-Search 是一个面向 RNA 同源序列搜索的项目。当前仓库包含三部分能力：

- RNA 语言模型推理与 contact map 预测：`RNA-ESM2/`
- contact map 构建、比对和序列搜索算法：`profile-contact-alignment/`
- 数据库存储与检索接口：规划中，当前数据库文件暂未放入仓库

典型流程是：先用 RNA-ESM2 从 FASTA 序列预测 profile/contact map，再用 `profile-contact-alignment` 中的 C++ 程序构建或比对 contact map，最后根据比对分数筛选同源候选序列。

## 目录结构

```text
RESM-Search/
+-- RNA-ESM2/                    # RNA-ESM2 模型、推理脚本、训练脚本和示例 FASTA
|   +-- rna_esm/                  # 本地 RNA-ESM/ESM 模型代码
|   +-- evo/                      # MSA、tokenization、metrics 等辅助模块
|   +-- fasta/                    # 示例/测试序列
|   +-- ckpt/                     # 本地 checkpoint，默认不提交 Git
|   +-- data/                     # 本地预训练权重，默认不提交 Git
+-- profile-contact-alignment/    # contact map 比对与搜索 C++ 工具
    +-- apps/                     # 命令行入口
    +-- src/                      # 核心算法实现
    +-- include/                  # 头文件
    +-- scripts/                  # 构建、测试和 benchmark 脚本
    +-- docs/                     # 子项目文档
```

## 环境要求

Python 侧：

- Python 3.8+
- PyTorch
- pytorch-lightning
- hydra-core / omegaconf
- biopython
- numpy
- scipy
- scikit-learn
- pandas
- tqdm
- numba
- transformers
- tape-proteins
- fair-esm 或本地可用的 `esm` 包

C++ 侧：

- C++11 编译器
- CMake 3.0+
- pkg-config
- ViennaRNA / RNAlib2
- OpenMP，推荐开启

如果 `pkg-config` 找不到 RNAlib2，需要设置 ViennaRNA 的 pkgconfig 路径：

```bash
export PKG_CONFIG_PATH=/path/to/ViennaRNA/lib/pkgconfig:$PKG_CONFIG_PATH
```

## 权重和数据准备

大模型权重、checkpoint、搜索数据库通常不适合直接提交 GitHub。本仓库的 `.gitignore` 默认排除这些文件：

- `RNA-ESM2/ckpt/*.pt`
- `RNA-ESM2/ckpt/*.ckpt`
- `RNA-ESM2/data/*.pt`
- `*.npz`, `*.npy`, `*_allmaps.bpmap`
- C++ 编译产物和 Python 缓存

运行模型推理前，请在本地准备所需权重。代码中常见默认路径包括：

```text
RNA-ESM2/ckpt/RUNGPURNA-ESM2-trans-mappro-KDNY-epoch-30-step-15650.pt
RNA-ESM2/ckpt/RNA-ESM2-trans-2a100-mappro-KDNY-epoch=07-valid_F1=0.561-v1.ckpt
RNA-ESM2/data/esm2_t30_150M_UR50D.pt
RNA-ESM2/data/esm2_t33_650M_UR50D.pt
```

数据库目录暂未纳入仓库。用户需要把待搜索的 genome/database FASTA 放到本地路径，并在命令中用 `-b`、`-genome` 或 `-database` 指定。

## RNA-ESM2 推理

进入模型目录：

```bash
cd RNA-ESM2
```

预测 contact map 的常用脚本有：

- `esm2embatt_650M_con.py`：650M 配置，输出 contact map
- `esm2embatt_150M_con.py`：150M 配置，输出 contact map
- `esm_profile_test_batch.py`：批量测试/导出 profile 和 contacts
- `esm_profile_test_cpu.py`：CPU 推理版本

Hydra 参数可以直接在命令行覆盖：

```bash
CUDA_VISIBLE_DEVICES=0 python esm2embatt_650M_con.py \
  produce.seqfile=./fasta/testRNA42.fasta \
  produce.outname=testRNA42 \
  data.device=cuda:0 \
  data.model_path=./ckpt/RUNGPURNA-ESM2-trans-mappro-KDNY-epoch-30-step-15650.pt
```

常见输出：

- `<outname>_allmaps.bpmap`：所有序列的 contact map 文本
- `<outname>_mapdict.json`：序列 ID 到 map 的索引信息，部分脚本生成
- `<outname>_seqdict.json`：序列 ID 到序列的映射，部分脚本生成
- `<outname>_contacts.npz`：contact tensor，部分脚本生成
- `<outname>_features.npz`：profile/feature tensor，部分脚本生成

注意：部分历史脚本内部带有硬编码输出目录，例如 `/mnt/remote_home/...`。换机器运行前，请检查脚本中的 `results_path`，或通过命令行参数/代码修改为本机可写目录。

## 模型训练

`RNA-ESM2/main.py` 是训练入口，使用 Hydra 管理参数。默认配置会读取：

- `data.ffindex_path=data/UniRef30_2020_02_a3m`
- `data.trrosetta_path=data/trrosetta`
- `data.trrosetta_train_split=train.txt`
- `data.trrosetta_valid_split=test.txt`

示例：

```bash
cd RNA-ESM2
python main.py \
  data.ffindex_path=/path/to/UniRef30_2020_02_a3m \
  data.trrosetta_path=/path/to/trrosetta \
  train.gpus=1 \
  train.max_epochs=100
```

训练数据体量通常较大，请放在本地数据目录或外部存储中，不建议提交到 GitHub。

## 构建 contact map 比对工具

进入 C++ 子项目：

```bash
cd profile-contact-alignment
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 8
```

二进制文件会生成到：

```text
profile-contact-alignment/bin/
```

主要命令：

- `cmapbuild`：从 DBN/结构文件构建 contact map
- `cmapsearch`：用 query 结构搜索 genome/database FASTA
- `RESM_Search`：query 结构和目标序列的搜索/比对入口
- `RESM_Search_align`：contact map/序列比对入口
- `RESM_Search_alignNCBI`：NCBI 场景入口
- `BPmap`：两个预计算 contact map 的直接比对

也可以使用子项目自带脚本：

```bash
scripts/build_and_test.sh
scripts/benchmark_performance.sh
python scripts/verify_optimizations.py
```

## 输入文件格式

DBN 文件使用三行一组：

```text
>sequence_id
AUGCUAGCUA
(((....)))
```

FASTA 文件使用标准 FASTA 格式：

```text
>sequence_id
AUGCUAGCUA
```

contact map 文件支持如下行：

```text
LEN 100
CON 3 45 1.0
PRF 0 A N 0.85 0.05 0.05 0.05
```

规则：

- `LEN` 或 `SIZE` 表示 map 长度
- `CON i j score` 表示第 `i` 和 `j` 个位置之间的接触
- `PRF` 是可选 profile 信息
- 如果没有 `PRF`，`BPmap` 会用 `N` 占位打印序列/profile 相关输出

## 使用示例

从 DBN 构建 cmap：

```bash
./bin/cmapbuild -dbn query.dbn -o query.cmp
```

从 DBN 和 MSA 构建 cmap：

```bash
./bin/cmapbuild -dbn query.dbn -a2m query.a2m -o query.cmp
./bin/cmapbuild -dbn query.dbn -sto query.sto -o query.cmp
```

搜索 genome FASTA：

```bash
./bin/cmapsearch -dbn query.dbn -genome genome.fa
```

搜索 database FASTA：

```bash
./bin/cmapsearch -dbn query.dbn -database database.fa
```

使用 query DBN 与 target FASTA 做 RESM 搜索：

```bash
./bin/RESM_Search -a query.dbn -b targets.fa -threads 8 -sep_cut 1
```

直接比对两个预计算 contact map：

```bash
./bin/BPmap -a query.bpmap -b target.bpmap -threads 8 -use_prf
```

常用参数：

- `-a`：query RNA 二级结构或 contact map 文件
- `-b`：target FASTA 或 contact map 文件
- `-threads` / `-num_threads`：OpenMP 线程数，默认 `8`
- `-gap_o`：gap opening penalty，默认 `-1`
- `-gap_e`：gap extension penalty，默认 `-0.01`
- `-sep_cut`：最小序列间隔过滤阈值，默认 `1`
- `-iter`：迭代次数，默认 `20`
- `-use_prf`：使用 profile 得分
- `-prf_w`：profile 得分权重
- `-use_gap_ss`：对二级结构区域使用 gap penalty 调整，默认开启
- `-gap_ss_w`：二级结构 gap penalty 权重
- `-silent`：减少日志输出

## 推荐工作流

1. 准备 query RNA 序列和结构，或准备 FASTA 序列。
2. 用 `RNA-ESM2` 预测 profile/contact map。
3. 用 `cmapbuild` 或预测得到的 `.bpmap` 生成 query map。
4. 准备待搜索的 genome/database FASTA。
5. 用 `cmapsearch`、`RESM_Search` 或 `BPmap` 进行比对。
6. 根据输出的 alignment score/BestScore 排序，筛选候选同源序列。

## 注意事项

- RNA 序列中的 `T` 会在部分读入逻辑中转换为 `U`。
- 部分脚本是研究阶段脚本，包含历史集群路径或 PBS 脚本；迁移环境时需要先检查输入、输出和 checkpoint 路径。
- 模型推理对显存要求较高。650M 模型建议使用 GPU；无 GPU 时使用 CPU 版本脚本或缩小 batch/序列长度。
- GitHub 不适合存放数 GB 模型权重和数据库。建议把权重放在对象存储、网盘、Hugging Face、GitHub Release 或机构内部存储，并在 README/配置中记录下载地址。

## 版本管理建议

建议提交到 GitHub 的内容：

- Python/C++ 源代码
- 构建脚本和轻量测试脚本
- 小型示例 FASTA/DBN
- README 和文档

不建议提交的内容：

- `.pt` / `.ckpt` 模型权重
- 大型 FASTA 数据库
- `.npz` / `.npy` 推理结果
- `__pycache__`
- CMake 构建目录、二进制文件、benchmark 产物
