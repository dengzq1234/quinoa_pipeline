# 藜麦耐盐性分析报告 — 回应 Raquel 的问题

日期：2026-09-01
分析流程：`quinoa_annotation_pipeline`（Step 1–5，全部 R 实现）
核心问题：**在 NaCl 胁迫下，哪些基因在耐盐（ST）与不耐盐（NT）品系之间的响应不同？这些基因的生物学功能是什么？**

---

## 0. 提纲：Raquel 的问题与回答路径

| # | Raquel 的问题 | 使用的分析 | 一句话结论 | 支撑文件 |
|---|---|---|---|---|
| **M** | rentrez + blastp 是否是做功能注释的正确方法？ | 用 biomaRt（Ensembl Plants）策展直系同源替代本地 BLAST | 存在更好的方案：biomaRt 一步映射到拟南芥直系同源，无需 FASTA / 本地 BLAST，全程 R | `pipeline_update_summary.txt` |
| **核心** | 哪些基因 ST vs NT 对 NaCl 响应不同？功能是什么？ | DESeq2 交互项 LRT（Step 1）+ GO 富集（Step 3/5） | 2,551 个差异基因，功能高度集中在**光合作用 / 叶绿体类囊体 / 色素（叶绿素）代谢 / 细胞壁**；ST 维持并上调光合机构，NT 则崩溃 | `degs_salt_tolerance/`, `05_cluster_go/` |
| **Q1** | 三个 qPCR 验证基因（PPR40 / AS-1 / trmH）是否受 NaCl 调控？ | 简化模型 `~ Tol + Treatment`，提取 NaCl 主效应（T30） | 三个基因**均无显著 NaCl 主效应**（padj 0.73 / 0.91 / 0.90），\|log2FC\| < 0.3 | `04_cluster/Q1_qPCR_NaCl_main_effect.tsv` |
| **Q2** | 差异响应基因能否按表达模式分组？ | `degPatterns`（DEGreport），top 1,000 DEG | 得到 **4 个聚类**，947 个基因；两大主轴 = ST 特异诱导（C2, n=495）与 NT 崩溃（C3, n=178） | `04_cluster/Q2_*` |
| **Q3** | 是否存在离群样本？ | 全样本 PCA，标记距组质心 > 2 SD 的样本（Tol × Treatment × Time 分组） | **未发现离群样本**（0 个） | `04_cluster/Q3_outlier_samples.tsv` |

---

## 1. 数据与方法

- **输入**：`gene_count.csv`（58,884 基因 × 73 样本，原始 read counts）+ `samples.csv`（`Treatment` = CONT/NaCl，`Tol` = ST/NT，`Time` = T0/T30）。
- **低表达过滤**：`rowSums(counts) ≥ 10` → 约 44,500 个表达基因。
- **分组变量**：使用 `Tol`（ST/NT 生物学分类），**不是** `Tolerance`（品系编号 14/227/287/460/577/670）。
- **差异表达模型（Step 1，仅 T30 胁迫后样本）**：
  - 全模型 `~ Tol + Treatment + Tol:Treatment`
  - LRT 对照缩减模型 `~ Tol + Treatment`
  - 提取交互项 `TolST.TreatmentNaCl`：正值 = NaCl 在 ST 中比在 NT 中更强诱导（或更少抑制）。
  - 显著性阈值：`padj < 0.05` 且 `|log2FC| ≥ 1`。
- **功能注释（Step 3）**：藜麦 NCBI gene ID → 藜麦 Ensembl ID → 拟南芥直系同源（biomaRt / Ensembl Plants，两段式查询避开 attribute-page 报错）。
- **GO 富集**：`clusterProfiler` + `org.At.tair.db`，BP / MF / CC 三本体，BH 校正。
- **聚类（Step 4）**：`degPatterns` 对 top 1,000 DEG（按 padj）做表达模式聚类，`time = Treatment`、`col = Tol`。
- **聚类 GO（Step 5）**：对每个聚类的拟南芥直系同源做 `enrichGO`，**以全部表达基因的直系同源集合（约 8,218 个拟南芥基因）为 universe**（`_corrected` 输出）。

---

## 2. 结果

### 2.1 方法学问题（M）：注释方案

原 Python/bash 流程用 `rentrez` 拉序列 + 本地 `blastp`/DIAMOND 打拟南芥。**已替换为 biomaRt（Ensembl Plants）策展直系同源映射**：

- 直接得到 quinoa → Arabidopsis 基因对，无需下载 FASTA、无需本地 BLAST；
- `clusterProfiler` 一步跑完 BP/MF/CC 三本体富集；
- FASTA 检索与 DIAMOND BLAST 保留为可选脚本（`scripts/optional/`），不再是主流程；
- 原 Python 脚本归档于 `scripts/legacy/`。

> 在 1,411 个未表征藜麦基因上的测试运行验证了端到端可用性（`gene_list_1411/`），富集到光合作用、叶绿体类囊体、糖苷水解酶活性等生物学上合理的条目。**注意该测试用的是背景基因表，不是 DEG 表**，只说明方法可行。

### 2.2 差异表达（Step 1）

| 指标 | 数值 |
|---|---|
| 显著 DEG（padj < 0.05，\|log2FC\| ≥ 1） | **2,551** |
| ST 中更被 NaCl 诱导（up in ST） | 2,337 |
| ST 中更被 NaCl 抑制（down in ST） | 214 |

输出：`01_deseq2/` 下的完整结果表、PCA、火山图、top 500 热图；DEG ID 列表 `03_DEG_list_for_annotation.txt` 作为 Step 3/4 的输入。
强烈的方向性不对称（2,337 vs 214）说明：**ST 的耐盐性主要来自"额外开启"一批基因，而非关闭基因。**

### 2.3 DEG 层面 GO 富集（Step 3）

2,551 个 DEG 中有 **637** 个映射到拟南芥直系同源（约 25%）。GO 富集（119 个显著条目，`degs_salt_tolerance/04_GO_enrichment_results.tsv`）的主线：

- **细胞壁 / 多糖 / 木聚糖**：`polysaccharide metabolic process`、`plant-type secondary cell wall biogenesis`、`xylan metabolic process`、`apoplast`（CC，p.adj ≈ 5×10⁻⁹）— CESA4/7、IRX9/15、XTH、GUX2、FLA11 等。
- **色素 / 叶绿素代谢**：`pigment metabolic process`（p.adj ≈ 2×10⁻⁵）、`pigment biosynthetic process`、`porphyrin` — CHLM、FLU、ACSF、LIL3、BPG1。
- **脂肪酸 / 蜡质 / 苯丙烷类**：`fatty acid biosynthetic process`（p.adj ≈ 2×10⁻⁶）、`phenylpropanoid biosynthetic process`、`very-long-chain fatty acid` — 角质层/表皮蜡相关（CER4/6/26、KCS、LACS2）。
- **叶绿体组织 / 蛋白定位到叶绿体**：`plastid organization`（p.adj ≈ 9×10⁻⁸）、`protein localization to chloroplast`。
- **胁迫响应**：`response to water deprivation`、`response to cold`、`response to UV`、`antioxidant activity`、`detoxification`（GST、过氧化物酶、SOD/FSD3）。

### 2.4 Q1 — qPCR 验证基因的 NaCl 主效应

模型 `~ Tol + Treatment`（T30），提取 `Treatment_NaCl_vs_CONT`：

| gene_id | 基因 | baseMean | log2FC | pvalue | padj | 显著? |
|---|---|---|---|---|---|---|
| 110697655 | PPR40 | 120.4 | −0.275 | 0.589 | 0.732 | 否 |
| 110710750 | AS-1 | 241.2 | +0.110 | 0.833 | 0.905 | 否 |
| 110703716 | trmH | 63.4 | +0.175 | 0.818 | 0.895 | 否 |

**结论**：在 T30、把 ST 与 NT 合并看的情况下，这三个基因对 NaCl **没有显著的整体转录响应**，效应量极小（\|log2FC\| < 0.3）。
提示：如果 qPCR 里看到了变化，可能是 (a) 品系特异（需要看交互项而非主效应）、(b) 时间点差异（T0 vs T30 或更早）、或 (c) 转录后调控。当前脚本只报告了 NaCl 主效应。

### 2.5 Q2 — 表达模式聚类

`degPatterns` 对 top 1,000 DEG 聚类，947 个基因落入 **4 个聚类**：

| 聚类 | 基因数 | 表达模式 |
|---|---|---|
| **1** | 257 | 交叉型：NaCl 下 ST↑ / NT↓ |
| **2** | 495 | ST 特异强诱导（NT 基本不动） |
| **3** | 178 | NaCl 下 NT 崩溃，ST 保持稳定 |
| **4** | 17 | ST 特异下调（NT 稳定） |

C2 + C3 合计占 71%，构成同一生物学轴的两极：**ST 主动上调 / NT 被动丧失。**
输出：`04_cluster/Q2_degPatterns_clusters.png`、`Q2_cluster_assignments.tsv`。

### 2.6 Q3 — 离群样本

全样本 blind VST + PCA，按 `Tol × Treatment × Time` 分组，标记距组质心 > 2 SD 的样本。
**结果：0 个离群样本**（`Q3_outlier_samples.tsv` 为空）。数据质量良好，无需剔除样本。
输出：`04_cluster/Q3_PCA_outliers.png`。

### 2.7 各聚类 GO 富集（Step 5，genome-wide universe）

以约 8,218 个"表达基因直系同源"为 universe（比只用 DEG 做背景更保守、更正确）：

| 聚类 | 显著条目数 | 顶端条目 | 生物学解读 |
|---|---|---|---|
| **1**（交叉 ST↑/NT↓） | 2 | `phototropism`（p.adj 0.019；RPT2/JK218/JK224）、`response to lipid`（p.adj 0.043；LOX2/LTP3/annexin/GASA） | ST 与 NT 在光信号与膜/脂质信号上分化 |
| **2**（ST 特异诱导） | 40 | **`thylakoid`（CC，p.adj 4×10⁻³³）**、`photosynthesis`（BP，p.adj 2×10⁻¹⁷）、`photosynthesis, light harvesting`、`photosynthetic membrane`、`plastid stroma/envelope` | **ST 在 NaCl 下主动上调叶绿体光捕获与光合机构**（LHCA1/2/3、LHCB2/3、CAB4、CP24、PSB*、PSA*、PETC、RBCS1A、SBPASE、PRK） |
| **3**（NT 崩溃） | 23 | `pigment metabolic process`（p.adj 3×10⁻⁵）、`tetrapyrrole / porphyrin metabolic process`、`chlorophyll metabolic process`、`lipid / wax biosynthetic process`、`cell wall`（CC） | **NT 在 NaCl 下丧失叶绿素/色素合成能力与细胞壁完整性**（AtCLA1、ACLA-1、GBP、PSB29、CER6、KCS、XTH9、GASA1、BGAL3） |
| **4**（ST 下调） | — | 直系同源基因 < 5，跳过富集 | 样本量不足，无法解读 |

---

## 3. 生物学结论

1. **耐盐性的核心是"保住光合作用"。** 差异全部围绕叶绿体：ST 品系在 NaCl 胁迫下*上调*类囊体光捕获复合体和光合电子传递链（Cluster 2），而 NT 品系*失去*叶绿素/四吡咯生物合成能力（Cluster 3）。二者是同一生理轴的正负两极。
2. **ST 靠"加法"而非"减法"。** 2,337 个基因在 ST 中更被诱导，仅 214 个更被抑制；表达模式上 ST 特异诱导（C2）远多于 ST 特异下调（C4，仅 17 个）。
3. **结构层面的第二条线索是细胞壁与角质层。** DEG 层面与 Cluster 3 都富集到次生细胞壁生物合成（CESA/IRX/XTH）和超长链脂肪酸/蜡质（CER/KCS）——NT 可能同时在细胞壁与表皮屏障上失守。
4. **三个 qPCR 基因（PPR40/AS-1/trmH）不是这条主线上的驱动基因**——至少在 T30 的 NaCl 主效应层面看不到响应。若要做验证实验，Cluster 2 的 LHC 基因和 Cluster 3 的叶绿素合成基因是更有说服力的候选。

---

## 4. 局限与后续

- **注释覆盖率**：仅约 19% 的表达藜麦基因、约 25% 的 DEG 能映射到拟南芥直系同源。GO 结果偏向"已知、保守"的通路，藜麦特异或未表征基因被系统性漏掉。
- **`_corrected` 聚类 GO 的可复现性**：已修复。`scripts/05_cluster_go_enrichment.R` 现在直接读取 Step 5a 的 `background_orthologs_all.tsv`（约 8,218 个拟南芥基因）作为 GO universe，输出 `*_GO_corrected.tsv` / `*_GO_corrected_dotplot.png` / `all_clusters_GO_corrected_combined.tsv`；缺背景文件会报错提示先跑 5a。旧的 `*_GO_results.tsv`（DEG-only 背景）保留为历史产物。
- **Q1 用的是 NaCl 主效应而非交互项**：如果 Raquel 关心的是"这三个基因在 ST/NT 之间响应是否不同"，应补充报告它们的 `TolST.TreatmentNaCl` 交互项结果。
- **流程整合**：Step 4、5、5a 目前独立运行，尚未接入 `run_pipeline.R`。
- **路径**：`run_pipeline.R` 与各脚本顶部的默认路径仍指向 `/home/ziqi/...`，换机器需修改。
- **README 目录名与实际不一致**：README 写 `results/01_deseq2/`、`results/03_annotation/`，实际 DEG 层面输出在 `results/degs_salt_tolerance/`。

### 建议的下一步
1. 补 qPCR 三基因的交互项结果 + T0/T30 时间轴对比图。
2. 对 Cluster 2/3 的 LHC 与叶绿素合成基因做靶向可视化（箱线图：ST/NT × CONT/NaCl），作为图 1 候选。
3. 考虑用藜麦自身基因组注释（若有）或 eggNOG-mapper 补充直系同源，评估覆盖率提升。
4. 把 Step 4/5/5a 接入 `run_pipeline.R`，统一路径为相对路径。

---

## 5. 文件索引

```
results/
├── degs_salt_tolerance/
│   ├── 04_quinoa_arabidopsis_orthologs_2551.tsv   # 637 对 quinoa→Arabidopsis
│   ├── 04_GO_enrichment_results.tsv               # DEG 层面 GO（119 条）
│   └── 04_GO_dotplot.png
├── 04_cluster/
│   ├── Q1_qPCR_NaCl_main_effect.tsv               # PPR40 / AS-1 / trmH
│   ├── Q2_degPatterns_clusters.png                # 4 聚类表达模式图
│   ├── Q2_cluster_assignments.tsv                 # 947 基因 × 聚类
│   ├── Q3_PCA_outliers.png
│   └── Q3_outlier_samples.tsv                     # 空 = 无离群样本
├── 05_cluster_go/
│   ├── background_orthologs_all.tsv               # GO universe（~8,218 拟南芥基因）
│   ├── cluster1_GO_corrected.tsv / _dotplot.png   # 2 条
│   ├── cluster2_GO_corrected.tsv / _dotplot.png   # 40 条（光合 / 类囊体）
│   ├── cluster3_GO_corrected.tsv / _dotplot.png   # 23 条（色素 / 细胞壁）
│   └── all_clusters_GO_corrected_combined.tsv     # 65 条合并
├── gene_list_1411/                                # 方法学测试运行（非 DEG）
└── pipeline_update_summary.txt                    # 2026-07-14 流程迁移说明
```
