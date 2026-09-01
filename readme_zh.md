# 藜麦耐盐分析流程

端到端 RNA-seq 分析流程，用于识别耐盐（ST）与非耐盐（NT）藜麦品系在
NaCl 胁迫下表达差异的生物学功能。

**核心生物学问题：** 哪些基因在 ST 和 NT 品系中对 NaCl 的响应不同？这些基因的功能是什么？

---

## 目录结构说明

本流程依托两个并列目录运行：

```
quinoa_raquel/
├── rnaseq_analysis/          # 上游：原始数据 + DESeq2 差异表达分析
│   ├── data/
│   │   ├── gene_count.csv    # 58,884 基因 × 73 样本（原始计数矩阵）
│   │   └── samples.csv       # 样本元数据
│   ├── scripts/
│   │   └── deseq2_analysis.R # 原始 DESeq2 脚本
│   └── results/
│       ├── 02_DESeq2_T30_TolxTreatment_full.tsv
│       └── 03_DEG_list_for_annotation.txt   # 2,551 个 DEG → 步骤 3 和 4 的输入
└── quinoa_pipeline/          # 本流程：注释 + 聚类 + GO 富集
    ├── scripts/
    └── results/
```

`quinoa_pipeline` 从 `rnaseq_analysis/` 读取原始数据和 DEG 列表，
所有输出写入自身的 `results/` 目录。

---

## 环境要求

首次运行前，安装所有必要的 R 包：

```r
source("scripts/00_install_packages.R")
```

---

## 项目结构（quinoa_pipeline/）

```
├── data/
│   ├── input/              # 独立基因 ID 列表（仅注释时使用）
│   └── reference/          # 拟南芥 DIAMOND 数据库（可选 BLAST 步骤）
├── results/
│   ├── degs_salt_tolerance/ # 同源注释 + GO 富集（2,551 DEG 级别）
│   ├── 04_cluster/          # 表达聚类（degPatterns）+ 离群样本 QC
│   └── 05_cluster_go/       # 每个表达簇的 GO 富集
├── scripts/
│   ├── 00_install_packages.R       # 一次性安装
│   ├── 01_deseq2_analysis.R        # 步骤 1：ST vs NT 差异表达
│   ├── 02_filter_gene_list.R       # 步骤 2（可选）：过滤未注释基因
│   ├── 03_annotate_and_go.R        # 步骤 3：biomaRt 同源映射 + GO 富集
│   ├── 04a_qpcr_gene_expression.R  # 步骤 4a：qPCR 基因 NaCl 响应检验
│   ├── 04_cluster_analysis.R       # 步骤 4：表达聚类 + 离群样本检测
│   ├── 05a_build_background_orthologs.R  # 步骤 5a：构建 GO 背景集
│   ├── 05_cluster_go_enrichment.R  # 步骤 5：每簇 GO 富集
│   ├── optional/
│   │   ├── get_fasta.R             # 从 NCBI 获取蛋白质 FASTA
│   │   └── blast_against_arabidopsis.sh  # DIAMOND BLAST vs 拟南芥
│   └── legacy/                     # 原始 Python 流程（归档）
└── run_pipeline.R                  # 步骤 1–3 入口
```

---

## 运行方式

编辑 `run_pipeline.R` 顶部的两个路径：

```r
rnaseq_data_dir <- "/path/to/rnaseq_analysis/data"   # gene_count.csv + samples.csv
results_base    <- "/path/to/quinoa_pipeline/results"
```

在 `quinoa_pipeline/` 目录下依次运行：

```r
source("run_pipeline.R")                              # 步骤 1–3
source("scripts/04a_qpcr_gene_expression.R")          # 步骤 4a：qPCR 基因检验
source("scripts/04_cluster_analysis.R")               # 步骤 4：聚类 + 离群 QC
source("scripts/05a_build_background_orthologs.R")    # 步骤 5a（仅需运行一次）
source("scripts/05_cluster_go_enrichment.R")          # 步骤 5：簇级 GO 富集
```

步骤 4、5a、5 均可独立运行——每个脚本都有 `if (!exists(...))` 路径守卫，
可单独 source 或从主脚本调用。

---

## 各步骤说明

### 步骤 1 — DESeq2 差异表达（`01_deseq2_analysis.R`）

**输入：** `gene_count.csv`（58,884 基因 × 73 样本）、`samples.csv`

**做什么：**
- 过滤到 T30 样本（盐胁迫后期）
- 模型：`~ Tol + Treatment + Tol:Treatment`，LRT 似然比检验
- 检验哪些基因在 ST 和 NT 品系中对 NaCl 的响应不同
- 使用 `Tol`（ST/NT 生物学分类），**不是** `Tolerance`（品系编号）

**输出（`rnaseq_analysis/results/`）：**
- `01_PCA_all_samples.png` — 全部 73 样本的 QC PCA 图
- `02_DESeq2_T30_TolxTreatment_full.tsv` — 含基因注释的完整结果表
- `03_DEG_list_for_annotation.txt` — DEG 基因 ID 列表（步骤 3 和 4 的输入）
- `03a_DEGs_up_in_ST.tsv` — 2,337 个 NaCl 下 ST 上调更多的基因
- `03b_DEGs_down_in_ST.tsv` — 214 个 NaCl 下 ST 下调更多的基因
- `04_PCA_T30.png` — T30 样本 PCA 图
- `05_volcano_T30_ST_vs_NT.png` — 火山图
- `06_heatmap_top500.png` — 前 500 个 DEG 热图

---

### 步骤 2 — 过滤未注释基因（`02_filter_gene_list.R`）— 可选

**输入：** 基因 ID 列表（每行一个 NCBI gene ID）

**做什么：** 通过 `rentrez` 查询 NCBI 获取基因元数据，保留未注释的蛋白编码基因。
默认在 `run_pipeline.R` 中跳过，因为 GO 富集可直接在完整 DEG 列表上运行。

> 注意：约需 17 分钟（2,500 个基因），受 NCBI API 速率限制（0.4 秒/基因）。

**输出：** TSV 文件，列：`gene_id, gene_name, gene_desc, gene_type,
refseq_genomic, refseq_mrna, refseq_peptide, related_acc`

---

### 步骤 3 — 拟南芥同源映射 + GO 富集（`03_annotate_and_go.R`）

**输入：** `rnaseq_analysis/results/03_DEG_list_for_annotation.txt`（2,551 个 DEG）

**做什么：**
- 藜麦 NCBI gene ID → 藜麦 Ensembl ID → 拟南芥同源基因
  （通过 Ensembl Plants biomaRt，两次查询避免属性页冲突错误）
- 使用 `clusterProfiler` + `org.At.tair.db` 跑 BP、MF、CC 三个本体的 GO 富集

**输出（`results/degs_salt_tolerance/`）：**
- `04_quinoa_arabidopsis_orthologs_2551.tsv` — 藜麦→拟南芥同源对（637 对）
- `04_GO_enrichment_results.tsv` — 完整 GO 富集结果表
- `04_GO_dotplot.png` — 每个本体 top 15 GO 条目气泡图

---

### 步骤 4a — qPCR 基因表达检验（`04a_qpcr_gene_expression.R`）

**输入：** `gene_count.csv`、`samples.csv`

**做什么：**
- 检验 Raquel 的三个 qPCR 验证基因在 T30 阶段是否对 NaCl 显著响应
- 模型：`~ Tol + Treatment`（NaCl 主效应）
- 检验基因：PPR40（110697655）、AS-1（110710750）、trmH（110703716）

**输出（`results/04_cluster/`）：**
- `Q1_qPCR_NaCl_main_effect.tsv` — 三个基因的 log2FC、p-value、padj

**结论：** 三个基因均无显著差异表达（padj > 0.7），在各条件下表达稳定，
适合作为 RT-qPCR 内参基因候选。

---

### 步骤 4 — 表达聚类 + 离群样本 QC（`04_cluster_analysis.R`）

**输入：** `gene_count.csv`、`samples.csv`、步骤 1 产生的 DEG 列表

**做什么：**
- **Q2** — 用 `degPatterns`（DEGreport）对 padj 最小的前 1,000 个 DEG 做表达聚类，
  识别在 ST/NT × NaCl 条件下具有共同表达轨迹的基因组
- **Q3** — 离群样本检测：对全部样本做 PCA，标记距组质心（Tol × Treatment × Time）
  超过 2 SD 的样本

**输出（`results/04_cluster/`）：**
- `Q2_degPatterns_clusters.png` — 表达轨迹图（共 4 个簇）
- `Q2_cluster_assignments.tsv` — 每基因的簇分配（947 基因 × 4 簇）
- `Q3_PCA_outliers.png` — 标注离群样本的 PCA 图
- `Q3_outlier_samples.tsv` — 被标记的离群样本列表

**簇汇总（共 947 个基因）：**

| 簇 | 基因数 | 表达模式 |
|----|--------|---------|
| 1 | 257 | 交叉型：NaCl 下 ST↑ / NT↓ |
| 2 | 495 | ST 特异强诱导：NaCl 下 ST 大幅上调，NT 无明显变化 |
| 3 | 178 | NT 崩溃型：NT 在 NaCl 下骤降，ST 稳定 |
| 4 | 17  | ST 特异下调：NaCl 下仅 ST 压制 |

---

### 步骤 5a — 构建 GO 背景集（`05a_build_background_orthologs.R`）

**在步骤 5 之前运行一次。** 输出有缓存，文件已存在则跳过。

**做什么：**
- 从 Ensembl Plants 下载完整藜麦基因表（biomaRt，不设过滤条件）
- 过滤到所有表达基因（rowSums ≥ 10，n = 44,524）
- 下载完整藜麦→拟南芥同源表（不设过滤条件）
- 本地 join，产生正确的 GO 富集背景集

> **技术说明：** biomaRt 2.66 的 httr2 POST 后端与 Ensembl Plants 不兼容
>（HTTP 405）。本脚本通过 `assignInNamespace` 将 `.submitQueryXML` 替换为
> GET 请求，并下载完整表后在本地过滤，绕过过滤查询失败的问题。

**输出：** `results/05_cluster_go/background_orthologs_all.tsv`  
8,401 条表达藜麦基因→拟南芥同源对 → 8,218 个唯一拟南芥基因作为 GO 富集背景集。

---

### 步骤 5 — 每簇 GO 富集（`05_cluster_go_enrichment.R`）

**输入：**
- 步骤 4 的簇分配（`Q2_cluster_assignments.tsv`）
- 步骤 3 的 DEG 同源表（`04_quinoa_arabidopsis_orthologs_2551.tsv`）
- 步骤 5a 的背景集（`background_orthologs_all.tsv`）

**做什么：**
- 对每个簇，提取对应的拟南芥同源基因
- 用全体表达基因同源集作为背景（而非仅 DEG），运行 `enrichGO`
  （clusterProfiler，org.At.tair.db，keyType="TAIR"，ont="ALL"）
- 保存每簇的 TSV 结果和气泡图

**输出（`results/05_cluster_go/`）：**
- `cluster1_GO_corrected.tsv` / `_dotplot.png` — 2 个显著条目
- `cluster2_GO_corrected.tsv` / `_dotplot.png` — 40 个显著条目
- `cluster3_GO_corrected.tsv` / `_dotplot.png` — 23 个显著条目
- `all_clusters_GO_corrected_combined.tsv` — 65 个条目合并表

**主要生物学发现：**

| 簇 | 顶级 GO 发现 | 生物学解读 |
|----|-------------|-----------|
| 2 | 类囊体膜（CC，p.adj = 4×10⁻³³）；光合作用（BP，p.adj = 2×10⁻¹⁷） | ST 在 NaCl 下主动上调叶绿体光捕获机器 |
| 3 | 色素代谢；卟啉/四吡咯代谢；细胞壁 | NT 在 NaCl 下丧失叶绿素合成能力和细胞壁完整性 |
| 1 | 向光性；脂质响应 | ST 和 NT 在光信号和膜信号上的分叉响应 |

簇 2 和簇 3 描述同一生物学轴的两个极端：ST 维护/诱导光合能力，
而 NT 同时失去色素合成和结构完整性。

---

## 可选步骤

仅在需要蛋白质 FASTA 序列时使用（如系统发育、结构预测）：

```r
source("scripts/optional/get_fasta.R")
# bash scripts/optional/blast_against_arabidopsis.sh
```

---

## 注意事项

- `Tol`（ST/NT）是正确的生物学分组变量，**不是** `Tolerance`
  （品系编号：14、227、287、460、577、670）
- GO 富集使用拟南芥注释，因为藜麦缺乏完整的 GO 数据库；
  约 19% 的表达藜麦基因可映射到拟南芥同源基因
- 步骤 4 和 5 独立于 `run_pipeline.R` 运行（尚未整合进主流程）
- 同源映射使用 Ensembl Plants 的精选映射，而非 BLAST 比对
- 步骤 5a 中的 biomaRt GET patch 仅适用于 biomaRt ≥ 2.66 在
  Ensembl Plants 上的使用；早期版本或其他 mart 不需要
- `data/input/gene_list_1411.txt` 是早期批次数据的候选基因列表，
  当前流程不再使用，仅作历史参考保留。当前 DEG 列表为 2,551 基因，
  来自 `rnaseq_analysis/results/03_DEG_list_for_annotation.txt`
