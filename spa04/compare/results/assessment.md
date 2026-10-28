# Space Ranger 结果比较

比较 `spa04/Visium_HD_Mouse_Brain`（Space Ranger spaceranger-4.0.1）与 `spa04/out_10x`（Space Ranger spaceranger-3.0.0）。
两份输出报告相同的 `334,184,277` 条 reads。

## 结论

- 两份结果都对应 Mouse Brain Visium HD CytAssist 样本和 slide `H1-8CPZTMY`；mapping/QC 指标非常接近。
- tissue mask Jaccard：2 µm `0.9873`，8 µm `0.9862`，16 µm `0.9865`。差异可见 `tissue_mask_comparison.png`。
- 8 µm raw bin UMI 总量相关 `0.999977`，全基因 pseudobulk UMI 相关 `0.999999943`；整体信号高度一致，但逐基因逐 bin 矩阵不是完全相同。
- 输出中的 CytAssist/组织高低分辨率图像 SHA256 核验：`3/3` 个文件逐字节相同；详见 `image_comparison.csv`。
- 本地 4.0.1 结果报告 `40223` 个细胞且有 segmentation 文件；官网 3.0.0 输出没有相同口径的细胞级结果。
- 最重要的流程差异：官网下载 HTML 命令包含 `--loupe-alignment`，本地 4.0.1 invocation 未使用手工 alignment；tissue mask 与 bin 数的差异不能直接解释为表达定量失败。

## 关键 QC

| 指标 | 本地 4.0.1 | 官网 3.0.0 | 本地减官网 |
|---|---:|---:|---:|
| Valid Barcodes | 91.9790% | 91.6220% | +0.3570 个百分点 |
| Sequencing Saturation | 70.0598% | 70.0732% | -0.0134 个百分点 |
| Reads Mapped Confidently to Probe Set | 98.7388% | 98.7388% | +0.0000 个百分点 |
| Fraction Reads in Squares Under Tissue | 98.4267% | 98.3214% | +0.1053 个百分点 |

### Raw bin UMI 与 tissue-bin 数

| bin | raw UMI 4.0.1 | raw UMI 3.0.0 | UMI 差异 | tissue bins 4.0.1 | tissue bins 3.0.0 |
|---|---:|---:|---:|---:|---:|
| 2 µm | 90,808,424 | 90,415,420 | +0.4347% | 6,336,843 | 6,296,688 |
| 8 µm | 90,808,424 | 90,415,420 | +0.4347% | 397,587 | 393,543 |
| 16 µm | 90,808,424 | 90,415,420 | +0.4347% | 99,948 | 98,917 |

## 差异原因与限制

- 软件版本不同（4.0.1 对 3.0.0）；官网命令使用 Mouse Transcriptome Probe Set v2.0。本地 `outs/probe_set.csv` 文件头同样标为 `Visium Mouse Transcriptome Probe Set v2.0`、mm10、2020-A，因此不是仅因文件名不同就能认定探针版本不同；官网原始探针文件未在项目内，无法证明探针序列逐字节相同。
- 3.0.0 使用 Loupe 手工对齐 JSON，本地 4.0.1 使用自动配准。不同 image registration/tissue detection 是组织 mask 和 tissue-bin 数差异的合理候选原因；未控制变量重跑，不能将各类差异分别定量归因。
- raw UMI 总量也有 +0.4347% 的差异，不应全部归因于组织检测；barcode/UMI 处理的软件版本变化也是候选原因，但本次未核实版本变更文档。
- 矩阵的 ordered feature ID/名称逐项比较；barcodes 按 array row/column 空间坐标对齐，缺失矩阵列按零计数。bin UMI 和 pseudobulk 相关见 `matrix_comparison.csv`；不声称全量 gene×bin 计数完全一致。
- 4.0.1 raw 矩阵保留完整空间网格，3.0.0 raw 只包含其中一部分 barcode 列。2 µm 多出的 `4,036,061` 列中，仅 `14,883` 列有非零 UMI，其余是零 UMI 列；不能把 raw 列数增加理解为多出相应数量的组织 bin 或细胞。
- 相同 read 总数和相近 QC 支持二者来自同一公开数据集，但没有原始 FASTQ/reference/probe-set 两边文件的 SHA256，不能声称输入文件字节完全相同。
- 本比较仅读取现有输出，不重跑 Space Ranger、不改动原始结果或 Conda 环境。

## 直接复现教程

- 用户提供的 10x 页面提供该样本处理后的官方下载结果；`out_10x/Visium_HD_Mouse_Brain_web_summary.html` 内嵌此数据集实际 `spaceranger count` 命令，可直接查看输入参数。
- Space Ranger 官方 `count` 文档解释命令参数。针对这个数据集的运行说明见 `../tutorial.md`；手工运行 3.0.0 命令还需要官方下载的 Loupe alignment JSON，单有输出文件并不能重建该输入。

## 结果文件

- `metrics_comparison.csv`：QC 指标逐项对比。
- `matrix_comparison.csv`：6 个 resolution/type 组合的矩阵维度、UMI 总数、barcode/gene 对齐统计。
- `tissue_mask_comparison.csv`、`tissue_mask_comparison.png`：tissue 区域交集和图示。
- `run_metadata.csv`：软件版本、原始命令、slide/layout 和告警。
