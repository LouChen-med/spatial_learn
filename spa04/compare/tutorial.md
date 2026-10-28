# 直接使用本数据集的 Space Ranger 运行说明

## 先回答“有没有直接用这个数据集的教程”

这次没有在线核实到独立的、一步一步讲解这个具体 Mouse Brain 样本的教程：搜索工具没有返回可用正文，直接访问 10x 网站得到 HTTP 429。
不能把其他 Mouse Brain/Visium 数据集教程当作本样本的教程。

不过，你已经下载的 `../out_10x/Visium_HD_Mouse_Brain_web_summary.html` 保存了**这个具体数据集实际使用的完整命令**。
下面的说明直接依据该命令，以及你的 `../Visium_HD_Mouse_Brain/_cmdline`、`_invocation` 和 HDF5 文件内部版本，而不是猜测官网参数。

数据页：<https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he>。

## 官网这份下载结果是怎样跑出来的

从官方 HTML 提取的命令如下，路径保留官方的占位写法，仅改变换行排版：

```bash
spaceranger count \
  --id=Visium_HD_Mouse_Brain \
  '--description=Gene expression library of Mouse Brain (Visium HD) using the Mouse Whole Transcriptome Probe Set' \
  --cytaimage=/path/to/images/Visium_HD_Mouse_Brain_cytaimage.tif \
  --image=/path/to/images/Visium_HD_Mouse_Brain_microscope.tif \
  --slide=H1-8CPZTMY \
  --area=A1 \
  --transcriptome=/path/to/refdata-gex-mm10-2020-A \
  --probe-set=/path/to/spaceranger3.0/spaceranger-3.0.0/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv \
  --fastqs=/path/to/fastqs/Visium_HD_Mouse_Brain_fastqs \
  --loupe-alignment=/path/to/loupe_alignment/Visium_HD_Mouse_Brain.json \
  --create-bam=false
```

该下载包的所有矩阵 HDF5 属性 `software_version` 都记录 **spaceranger-3.0.0**。
这是一份运行参数记录，并不意味着以上占位路径在你的机器存在；不要直接复制执行。

## 你的运行与官方有什么不同

| 项目 | 你的运行 | `out_10x` 官方运行 |
|---|---|---|
| Space Ranger | 4.0.1 | 3.0.0 |
| FASTQ read 总数 | 334184277 | 334184277 |
| 参考目录名 | refdata-gex-mm10-2020-A | refdata-gex-mm10-2020-A |
| slide / area | H1-8CPZTMY / A1 | H1-8CPZTMY / A1 |
| 图像配准输入 | 未提供 Loupe JSON，`loupe_alignment_file = null` | 明确传入 `--loupe-alignment=...Mouse_Brain.json` |
| probe-set 文件名 | Visium_HD_Mouse_Brain_probe_set.csv | Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv |
| probe-set 版本标识 | 输出 `outs/probe_set.csv` 文件头确认 v2.0 / mm10 / 2020-A | HTML 命令指定 v2.0 / mm10-2020-A |
| BAM | 不生成 | 不生成 |
| 细胞分割 | 有 `segmented_outputs` 和 40223 个细胞 | 这份结果未提供细胞级输出 |

本地 `outs/probe_set.csv` 文件头确认同样是 **Visium Mouse Transcriptome Probe Set v2.0**，并不是只知道文件名而猜版本。
参考目录和探针版本标签一致，仍不能证明两边文件内容完全相同；需要对应输入文件的内容/哈希才能验证。
图像参数的角色不变：`--cytaimage` 是 CytAssist 图，`--image` 是显微镜 H&E 图，不能仅按自定义文件名猜测并交换。
原 Notebook 的参数解释表把这两个角色描述颠倒；本说明已经更正，未改动原 Notebook。

你两套输出中的 CytAssist、tissue hires/lowres 图像 SHA256 全部相同，所以目前没有“使用了另一张组织图片”的证据。
组织检测/对齐的输出 JPG 不同并不等于原始组织图不同。

## 如果想严格复现 `out_10x`

1. 准备相同 FASTQ、相同 mm10-2020-A reference、相同 Mouse Probe Set v2.0、相同 CytAssist/H&E 图像。
2. 获取官方用于这次运行的 `Visium_HD_Mouse_Brain.json`，它是 **Loupe alignment 输入**，不是任意文件名相似的 JSON；当前 `out_10x` 没有该文件。
3. 用与下载输出相同的 **Space Ranger 3.0.0** 以及上述命令。即使用同一版本，没有官方 alignment 文件也不构成同参数复现。
4. 使用新的、尚不存在的项目内输出目录；不要覆盖你已经跑好的 `Visium_HD_Mouse_Brain`。
5. 按 `results/assessment.md` 的标准检查 reads、mapping、barcode 交集、tissue mask 和矩阵计数，而不只比较 HTML 文件大小。

项目内现有 `out_10x` 没有 `Visium_HD_Mouse_Brain.json`；缺少它时不要杜撰 `--loupe-alignment` 路径，也不要用 `scalefactors_json.json` 代替。

如果继续使用 4.0.1，则应将它视为**新版本重新处理同一数据集**，不是对官网 3.0.0 的逐字节复现。
资源参数 `--localcores=20`、`--localmem=200`、输出目录以及是否使用 `sudo` 不是这两份结果生物学差异的充分解释；通常不需要为了跑分析而使用 `sudo`，你当前输出的确由 root 持有，但本次未修改其权限。

## 本次没有执行的操作

- 没有安装/切换 Space Ranger，没有下载新的大数据。
- 没有读取 `/data/...` 的项目外输入；仅使用项目内已有运行记录。
- 没有重跑 `spaceranger count`，没有改变 Conda、Git 或原始输出权限。
