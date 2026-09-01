# 向 Raquel 汇报 — 藜麦耐盐性分析

汇报日期：2026-09-02

分析背景：T30（胁迫后）样本，DESeq2 交互项模型 `~ Tol + Treatment + Tol:Treatment`（LRT，
缩减模型 `~ Tol + Treatment`），提取交互项 `TolST:TreatmentNaCl`。
"差异表达基因（DEG）" = 交互项 `padj < 0.05` 且 `|log2FC| ≥ 1` → **2,551 个基因**
（2,337 个在 ST 中更被 NaCl 诱导，214 个更被抑制）。

---

## 第一部分：回答 Raquel 邮件的三个问题

### Q1 — PPR40 / AS-1 / trmH 是否为差异表达基因？

三个基因：

| qPCR 名 | Phytozome | NCBI Gene ID | 有 qPCR 数据的条件 |
|---|---|---|---|
| PPR40 | AUR62042953 | 110697655 | C / NaCl |
| AS-1  | AUR62030471 | 110710750 | C |
| trmH  | AUR62044515 | 110703716 | C / NaCl |

**目前能确定的（已在本机核对）：**

- **三个基因都不在 top-1,000 交互 DEG（按 padj 排序）中** → 不是强差异响应基因。
- **NaCl 主效应（T30，模型 `~ Tol + Treatment`）三个都不显著**：

  | 基因 | log2FC | pvalue | padj |
  |---|---|---|---|
  | PPR40 | −0.275 | 0.589 | 0.732 |
  | AS-1  | +0.110 | 0.833 | 0.905 |
  | trmH  | +0.175 | 0.818 | 0.895 |

  效应量极小（|log2FC| < 0.3）。
- 三个基因**都没有策展的拟南芥直系同源**，所以不会出现在 GO / ortholog 结果里（与 DEG 与否无关）。

**尚待确认（需要完整 DESeq2 结果表）：**

- 是否为"弱 / 边缘 DEG"（交互项 padj < 0.05 但按 padj 排名在 1,000 ~ 2,551 之间）。
- 需要 `02_DESeq2_T30_TolxTreatment_full.tsv`（当前不在本机），查这三个 gene ID 的
  baseMean、交互项 log2FC、padj，才能给出确定答案。

**初步结论**：三个基因**极可能不是差异响应基因** —— 无显著 NaCl 主效应，也不在
强 DEG 中。这与它们在 qPCR 设计中作为相对稳定 / 弱响应基因的角色一致。

**建议补充分析（对齐 qPCR）**：qPCR 数据是 C 与 NaCl（PPR40、trmH）/ 仅 C（AS-1）。
可额外报告：ST 内 NaCl 效应、NT 内 NaCl 效应、CONT 条件下 ST vs NT，
并画标准化 count 箱线图（ST/NT × CONT/NaCl，如需要再加 T0/T30）。

### Q2 — 按表达谱聚类，看每个 cluster 在不同条件和时间点的表达变化

**是的，我们清楚这是哪种分析**：按基因表达轨迹聚类，再看每个 cluster 的平均表达如何
随条件、时间点变化。标准做法是 `degPatterns`（DEGreport 包）；同类方法还有
Mfuzz 软聚类、WGCNA 共表达模块。

**我们已经做了一个 T30 版本**（Step 4）：对 top 1,000 个交互 DEG 聚类，947 个基因 →
**4 个 cluster**：

| Cluster | n | 表达模式（仅 T30） |
|---|---|---|
| 1 | 257 | 交叉型：NaCl 下 ST↑ / NT↓ |
| 2 | 495 | ST 特异强诱导 |
| 3 | 178 | NaCl 下 NT 崩溃，ST 稳定 |
| 4 | 17  | ST 特异下调 |

**当前局限**：x 轴只有 Treatment（CONT→NaCl），**没有纳入 T0**。

**下一步（正是 Raquel 想要的）**：用全部时间点样本重建 VST 矩阵，`degPatterns` 中
`time = "Time"`、`col = "Treatment"`、按 `Tol` 分面（或用合并因子 `Treatment:Time` 作 x 轴），
让每个 cluster 的平均表达画成跨 **CONT-T0 / CONT-T30 / NaCl-T0 / NaCl-T30** 的轨迹，
并按 ST / NT 分开。

> 注意：T0 是施加 NaCl 之前，设计上 CONT-T0 ≈ NaCl-T0；轨迹主要体现"T0 基线 → T30 分化"。
> 这仍然有价值 —— 能区分"胁迫诱导的差异"与"本来就存在的基线差异"。

需要：原始 count 矩阵（当前不在本机）。

### Q3 — 部分 accession 有一个重复聚类不好，会影响审稿吗？（Pablo 也发现了）

**当前分析的问题**：我们自动化的离群检测按 `Tol × Treatment × Time` 分组（>2 SD
判定），报告 **0 个离群样本**。但这个分组把每个 ST / NT 下的 3 个 accession 合并了，
单个 accession 层面偏离的重复会被掩盖。Raquel 观察的是 **accession 层面**
（6 个：14 / 227 / 287 / 460 / 577 / 670）的聚类。

**下一步**：按 `accession × Treatment × Time` 重做 QC —— 样本-样本相关性热图 +
层次聚类 + 每组组内距离，定位到具体是哪个 accession 的哪个重复。

**关于"是否影响审稿"（给 Raquel 的判断要点）：**

1. 每组一个轻度偏离的重复在 RNA-seq 里**很常见，通常不致命**。
2. n = 3 / 组时，一个偏重复主要是**降低检验力**（抬高离散度估计），对 fold-change
   偏倚不大；DESeq2 的离散度收缩 + Cook's distance + independent filtering 能缓解极端情况。
3. 稳妥做法（建议直接写进稿子）：
   - **公开展示 QC**：PCA + 样本相关性热图作为补充图；
   - **量化**：是真"离群"还是只是"更吵"（组内距离、silhouette）；
   - **敏感性分析**：去掉这些重复重跑 DESeq2，展示关键结果稳定
     （DEG 重叠率高、光合 GO 依旧富集）—— 这是最有力的审稿答复；
   - 可选：RUVSeq / sva 建模隐藏变异，或把重复作为协变量 —— 仅在必要时；
   - **不要悄悄删样本**，记录决策依据。
4. 既然 Pablo 独立也看到了，建议**主动**在稿子里放一段 QC 说明 + 补充图，
   而不是等审稿人问。

---

## 第二部分：聚类的后续功能分析（已完成）

对第一部分 Q2 的 4 个 cluster，取每个 cluster 基因的拟南芥直系同源做 GO 富集
（`enrichGO`，BP / MF / CC）。**背景（universe）用"全部表达基因的直系同源集合"
（约 8,218 个拟南芥基因），不是只用 DEG** —— 这是统计上正确的背景。

### Cluster 2（ST 特异诱导，n=495）— 40 个显著条目

| 类型 | 顶端条目 | p.adjust |
|---|---|---|
| CC | **thylakoid（类囊体）** | 4×10⁻³³ |
| CC | plastid / chloroplast thylakoid membrane、photosynthetic membrane | 3×10⁻³³ ~ 6×10⁻³² |
| BP | **photosynthesis（光合作用）** | 2×10⁻¹⁷ |
| BP | photosynthesis, light harvesting（光捕获）、light reaction | 3×10⁻¹⁵ ~ 9×10⁻¹⁹ |
| CC | plastid stroma、chloroplast envelope | 9×10⁻²³ ~ 9×10⁻¹⁸ |

代表基因：LHCA1/2/3、LHCB2/3、CAB4、CP24、PSB\*、PSA\*、PETC、RBCS1A、SBPASE、PRK

→ **ST 在 NaCl 胁迫下主动上调叶绿体光捕获复合体和光合电子传递链。**

### Cluster 3（NT 崩溃，n=178）— 23 个显著条目

| 类型 | 顶端条目 | p.adjust |
|---|---|---|
| BP | **pigment metabolic process（色素代谢）** | 3×10⁻⁵ |
| BP | tetrapyrrole / porphyrin metabolic process（四吡咯 / 卟啉，叶绿素前体） | 7×10⁻⁵ |
| BP | **chlorophyll metabolic process（叶绿素代谢）** | 6×10⁻⁴ |
| BP | lipid / wax / very-long-chain fatty acid biosynthesis（脂质 / 蜡质） | 6×10⁻⁵ ~ 3×10⁻³ |
| CC | **cell wall（细胞壁）** | 1×10⁻³ |

代表基因：AtCLA1、ACLA-1、GBP、PSB29（色素 / 叶绿素）；AtCER6、ATKAS2、PAS2（蜡质 / 脂肪酸）；
XTH9、GASA1、BGAL3、ATGER3（细胞壁）

→ **NT 在 NaCl 下丧失叶绿素 / 色素合成能力，同时细胞壁与表皮蜡质屏障受损。**

### Cluster 1（交叉 ST↑/NT↓，n=257）— 2 个显著条目

- `phototropism`（向光性，p.adjust 0.019；RPT2 / JK218 / JK224）
- `response to lipid`（对脂质的响应，p.adjust 0.043；ATLOX2 / LTP3 / annexin / GASA1 / GASA4）

→ ST 与 NT 在光信号与膜 / 脂质信号上出现分化（信号偏弱，供参考）。

### Cluster 4（ST 下调，n=17）

直系同源基因不足 5 个，无法做 GO 富集。

---

## 第三部分：整体生物学结论

1. **耐盐性的核心是"保住光合作用"。** ST 上调类囊体光捕获与光合机构（Cluster 2），
   NT 失去叶绿素 / 四吡咯合成能力（Cluster 3）—— 同一生理轴的正负两极。
   这与 Raquel 观察到的"光合相关基因结果很合理"一致。
2. **ST 靠"加法"而非"减法"。** 2,337 个基因在 ST 中更被诱导 vs 仅 214 个更被抑制；
   ST 特异诱导聚类（495）远大于 ST 特异下调聚类（17）。
3. **第二条线索：结构屏障。** DEG 层面和 Cluster 3 都富集次生细胞壁生物合成
   （CESA / IRX / XTH）与超长链脂肪酸 / 蜡质（CER / KCS）—— NT 可能在细胞壁与
   表皮屏障上同时失守。
4. **三个 qPCR 基因（PPR40 / AS-1 / trmH）不在这条主线上**（见 Q1）。

---

## 第四部分：待办 / 需要的数据

**需要把以下文件拷到本机（或在有数据的机器上跑）才能完成：**

| 用途 | 需要的文件 |
|---|---|
| Q1 确定答案（查三个基因的交互项 padj） | `rnaseq_analysis/results/02_DESeq2_T30_TolxTreatment_full.tsv` 或 `03_DEG_list_for_annotation.txt` |
| Q2 时间点扩展聚类 | `rnaseq_analysis/data/gene_count.csv` + `samples.csv` |
| Q3 accession 层面 QC + 敏感性分析 | 同上 |

**拿到数据后的具体动作：**

1. **Q1**：grep 三个 gene ID → 报告 baseMean / 交互 log2FC / padj / DEG（是 / 否）；
   补 ST 内、NT 内 NaCl 效应 + 标准化 count 箱线图。
2. **Q2**：`degPatterns` 纳入 T0 + T30，每个 cluster 平均表达画跨
   CONT-T0 / CONT-T30 / NaCl-T0 / NaCl-T30 的轨迹，按 ST / NT 分面。
3. **Q3**：样本相关性热图 + 层次聚类（按 accession × Treatment × Time）定位偏离重复；
   去掉后重跑 DESeq2 做敏感性分析（DEG 重叠率、光合 GO 稳定性）；
   为稿子准备一段 QC 说明 + 补充图。
4. **Cluster 2 / 3 靶向图**：LHC 基因、叶绿素合成基因画 ST/NT × CONT/NaCl 箱线图作主图候选。

---

## 附：关键文件

| 内容 | 路径 |
|---|---|
| Q1 NaCl 主效应结果 | `results/04_cluster/Q1_qPCR_NaCl_main_effect.tsv` |
| Q2 聚类图 / 分配（T30 版） | `results/04_cluster/Q2_degPatterns_clusters.png`、`Q2_cluster_assignments.tsv` |
| Q3 离群样本 PCA（Tol 分组版） | `results/04_cluster/Q3_PCA_outliers.png` |
| 聚类 GO（各聚类 + 合并） | `results/05_cluster_go/cluster{1,2,3}_GO_corrected.tsv`、`all_clusters_GO_corrected_combined.tsv` |
| 聚类 GO 点图 | `results/05_cluster_go/cluster{1,2,3}_GO_corrected_dotplot.png` |
| 完整分析报告 | `results/analysis_report.md` |
