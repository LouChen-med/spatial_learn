# results

- `assessment.md`：主要发现、QC 表和分析限制；建议先阅读。
- `metrics_comparison.csv`：每个 QC 指标原值、绝对差、相对差；百分比指标的绝对差是百分点；旧版 `Valid UMIs` 对齐新版 `Valid UMI Sequences`。
- `matrix_comparison.csv`：三个分辨率、raw/filtered 共六组矩阵的 feature 一致性、barcode 交并集、UMI 总量、bin UMI Pearson 和 gene pseudobulk Pearson。
- `matrix_inventory.csv`：全部 12 个矩阵的形状、非零数、UMI 总量和 tissue 位置数。
- `run_metadata.csv`：软件版本、样本、slide、配准矩阵、HTML 原始命令、告警及细胞级输出文件列表。
- `image_comparison.csv`：输出 CytAssist、tissue hires/lowres 图像的 SHA256，检查两边文件是否字节相同。
- `tissue_mask_comparison.csv`：不同分辨率组织判定的交集、各自独有 bin 和 Jaccard。
- `tissue_mask_comparison.png`：灰色为两者均判为组织，橙色仅本地判为组织，蓝色仅官网判为组织，白色均不是。
- `validation.txt`：比较运行命令、已完成的验证与未验证的边界；没有环境/Git 修改。

CSV 使用 UTF-8 BOM，方便中文软件识别。结果由 `../scripts/compare_spaceranger_results.py` 生成，不含原始 FASTQ 或大表达矩阵副本。
`fontlist-v*.json` 是 Matplotlib 自动生成的字体缓存，不属于比较结果，可忽略。
