# Step 1 —— 真变异 vs 测序/mapping 错误 的分子连锁置信度模型

用 UMI/co-barcode **分子连锁特征**,区分一个 SNP 是**真（低频）变异 / 测序错误 / mapping 错误**,
在 HG002 GIAB 真值上训练与验证。这是 cLFR/cWGS "分子连锁变异 re-scorer" 的 compute-light 核心
（CPU、LightGBM）。完整论证见 `../memory/plan.md`。

---

## 为什么用分子特征（差异化所在）

标准 pileup caller 只看 read 计数,它**看不到** alt 支持是否**在独立分子内一致**。规律：
- **真变异**：被**多个独立分子（UMI/barcode family）**支持,且每个分子内部一致。
- **测序错误**：零星散落在单个 read 里,分子内不一致。
- **mapping 错误**：骑在低 MAPQ / soft-clip / repeat 邻近的 read 上。

`02_extract_features.py` 抽的正是这些**分子分辨**的信号（独立分子数、within-molecule 一致性、
MAPQ/BQ/strand/softclip/read 位置偏置）—— 这是标准 pileup caller 用不到的平台信号。

---

## 真值：Step 1 vs Step 2（决定 controls 够不够）

- **Step 1（本代码）**：纯 HG002 一个样本给两样金真值——
  - confident BED 内的 **hom-ref 位点却有 alt read = 错误**（负类）；
  - **GIAB SNV 位点 = 真变异**（正类,germline AF ~50%/100%）。
  - 足以做**误差抑制 + 校准置信度**。**不需要 titration。**
- **Step 2（以后）**：低 AF 的真变异需要 GIAB **titration**（HG002 + HG003/HG004 按已知比例掺）。
  纯 HG002 没有低 AF 正类，**别用 Step 1 的数据 claim 低频 sensitivity**。

---

## 流水线

```bash
# 改 run.sh 里的路径,然后：
bash run.sh
```

1. `01_make_candidates.py` —— confident BED 内扫 pileup → 打标签位点
   （真=GIAB SNV；错误=confident hom-ref 却有 alt read）。
2. `02_extract_features.py` —— 按 UMI/barcode 分子分组,逐位点抽分子连锁特征。
   分子 id 来源：`--molecule-source readname_regex`（cLFR read name 里的 barcode）
   或 `--molecule-source tag --molecule-tag BX`。
3. `03_train_eval.py` —— LightGBM 真-vs-错误,**按染色体留出测试集**,校准只在 train 切片上拟合,
   输出 PR/ROC/Brier + 按 VAF 分层的 PR-AUC + feature importance。

---

## Money chart（消融）

`run.sh` 训两次：`--feature-set all` vs `--feature-set baseline`
（baseline = dp/alt_reads/ref_reads/vaf,即 pileup caller 本来就有的）。
**PR-AUC 之差**（在**较高 VAF 的错误 = mapping artifact** 处最大,因为 VAF 本身分不开它们）
= 分子/mapping 特征带来的增量。这就是"分子连锁特征打过标准 pileup caller"的那张图。

> 合成数据验证：构造"高 VAF 的 mapping 错误"时,baseline PR-AUC=0.687,
> 加分子/MAPQ 特征后=1.000（delta +0.313）。真实数字要在服务器上跑 HG002 得出。

---

## 诚实护栏（已内置）

- **先按染色体 split,再做任何别的**；校准只在 train 上拟合,绝不碰 test。
- **split 之前不做 scale/SMOTE**（这正是旧 `sqanti3.py` 虚高指标的泄漏根源）。
- 纯 HG002 上 真=高 VAF、错误=低 VAF,单靠 VAF 就能分开一大半；所以**报消融 + VAF 分层,
  而不是一个头条 AUC**。
- Step 1 只 claim 误差抑制 + 校准。低频 sensitivity 归 Step 2。

---

## 第 04 步：如何通过 re-score 得到"SNP 正确的 isoform"

置信度模型不只是打分,它能直接产出**校正后的 isoform 序列**。思路 = **拿 re-scorer 当误差模型去 polish consensus**：

**流程**

1. cLFR consensus 流程照常产出 per-molecule consensus FASTA + 比对到 reference/transcriptome。
2. **取差异**：consensus 与 reference 不一致的每个位点 = 候选 SNP。
3. 对每个候选,用 `02_extract_features.py` 在**支持该 consensus 的 reads** 上抽同一套分子特征。
4. `04_apply_rescore.py` 加载训好的模型,逐 SNP 判定：
   - `KEEP`  `p_true ≥ 阈值` → 真转录本 SNP,**保留**；
   - `REVERT` `p_true < 阈值` → 测序/mapping 错误,**回退成参考碱基（纠正掉）**；
   - `EDIT`  A>G / T>C 且落在已知 RNA 编辑位点（`--editing-bed` REDIportal）→ **RNA 编辑,是真生物学,
     单独标注、绝不当错误回退**。
5. 用 KEEP 集生成校正后的 isoform FASTA：

```bash
python 04_apply_rescore.py --features consensus_features.tsv --model out/claimA_all/model.txt \
  --out-prefix out/corrected --threshold 0.5 --editing-bed REDIportal.bed

bcftools sort out/corrected.pass.vcf -Oz -o out/corrected.pass.vcf.gz
tabix -p vcf out/corrected.pass.vcf.gz
bcftools consensus -f REF.fa out/corrected.pass.vcf.gz > corrected_isoforms.fa
```

**这一步把 flagship 从"我给 SNP 打分"升级成"我产出 SNP 校正过的 isoform 序列"** —— 一个可交付的产品输出。

**ERCC 上的 money chart（校正价值）**：ERCC consensus 本应与已知 ERCC 序列逐碱基一致,任何残留 SNP
都是 consensus 错误。测**校正前 vs 校正后的 per-base 错误率 / 假 SNP 率**,下降幅度 = re-score 校正的价值。

**诚实 caveat**

- 阈值是 sensitivity↔specificity 的权衡：回退太狠会抹掉真实低支持 SNP（ASE 偏斜的位点）。
  用校准后的概率设 operating point,按用途调。
- **RNA 编辑必须显式处理**,否则会把真生物学"纠正"掉——所以有 `EDIT` 这一类。
- Step 1 模型是 **DNA（HG002）训练**的；用到 RNA/isoform 上时,应在 **ERCC（误差真值）上重新校准/微调**,
  方法可迁移但要在 RNA 数据上校准。

---

## 依赖

`pysam, pandas, numpy, scikit-learn, lightgbm`（见 `environment.yml`）；
`bcftools`（校正 FASTA 用）。真值 VCF 需 bgzip+tabix；参考需 `.fai`。