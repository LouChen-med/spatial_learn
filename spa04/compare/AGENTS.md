# Visium HD Mouse Brain 结果比较

- 项目根目录仍为 `/home/lcpyspa/spa`；本任务新增代码和结果限定在 `spa04/compare/`。
- 只读输入：`../Visium_HD_Mouse_Brain/outs/` 为用户 Space Ranger 4.0.1 运行；`../out_10x/` 为官网下载的 Space Ranger 3.0.0 结果，以 HDF5 内部版本和 HTML 命令为准。
- `scripts/` 放比较代码，`results/` 放输出，`tutorial.md` 说明这个具体数据集的复现参数与缺失输入。
- 使用既有 Conda 环境 `cosmx-spatial`，通过 `/home/lcpyspa/miniconda3/bin/conda run -n cosmx-spatial` 运行；不安装依赖、不修改 Conda/Git、不执行 `spaceranger count`。
- 原始 FASTQ、图像和大矩阵不复制、不修改；比较时先对齐 feature IDs 和空间 barcode，不依赖 filtered HDF5 的随机列顺序。
- 本次比较覆盖全部 2/8/16 µm raw/filtered HDF5，计算总 UMI、逐 bin UMI、基因 pseudobulk 和 tissue mask；不生成新聚类/注释，不把同一 bin 当作同一个细胞。
- 2026-10-28：核实本地 4.0.1 使用自动图像配准，官网下载 3.0.0 命令提供 `--loupe-alignment`；两边 read 总量均为 334184277。官网请求返回 HTTP 429，未能在线核实独立的同数据集教程；以本地官方 HTML 提取的运行命令为教学依据。
- 2026-10-28：完整比较 12 个矩阵并验证统计一致性；8 µm raw bin UMI Pearson 为 0.9999766，tissue-mask Jaccard 为 0.9861868，raw 总 UMI 增加约 0.4347%。三个输出组织图像 SHA256 均一致，本地输出 probe-set 文件头确认 Mouse v2.0 / mm10 / 2020-A。报告、图片、可复算代码和验证记录均保存在本目录。
- 2026-10-28：按用户本次明确授权，只读 `/data/04_ST_introduction/data/spaceranger/Visium_HD_Mouse_Brain_alignment_file.json`；不允许向 `/data` 写入。解码确认 393543 个 tissue 8 µm bins，与官网 2/8 µm tissue mask 逐项一致；A1 layout 和阵列到 CytAssist 变换也一致。文件没有 `cytAssistInfo`，应区分阵列定位/组织选择与显微镜图像手工配准。解读保存在 `results/alignment_file_interpretation.md`，教程同步更新；不重跑、不改 Conda/Git。
