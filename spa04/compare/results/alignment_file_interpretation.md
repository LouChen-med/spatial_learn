# Visium_HD_Mouse_Brain_alignment_file.json 解读

## 检查范围

- 日期：2026-10-28。
- 原文件：`/data/04_ST_introduction/data/spaceranger/Visium_HD_Mouse_Brain_alignment_file.json`。
- `/data` 仅按用户本次授权只读访问；未修改源文件、权限或目录，未执行 Space Ranger。
- 文件大小：118332 字节。
- 整个 JSON 文件的 MD5：`342b70135ed2a3376fd726191f4407d2`。
- 整个 JSON 文件的 SHA256：`9fcd80358685539d613e58258b0b2eb922fcfb1dff5261166233e9accb70ab3e`。

## 核心结论

这是 Visium HD 的 Loupe 对齐/组织选择文件，不是表达矩阵，也不是普通的图像缩放参数文件。
它保存了芯片编号与区域、芯片布局标定、阵列到图像的空间变换，以及哪些空间 bins 算作组织。

它的内容与当前 `spa04/out_10x` 官方 3.0.0 输出高度对应，而且有直接数值证据：

- 8 µm 组织 mask 的 702244 个 bin 逐项比较，差异为 0。
- 将此 mask 映射到 2 µm 网格后，11222500 个 square 的组织标签也与官网输出逐项一致，差异为 0。
- 解码的 A1 slide layout 与官网输出 layout 最大绝对差约 1.11e-16。
- 合成得到的 2 µm spot 中心到 CytAssist 像素坐标变换与官网输出最大绝对差约 1.14e-13，属于浮点误差。

因此它不是仅仅“文件名相似”的 JSON，可以作为本样本 `--loupe-alignment` 输入。
以上证明组织选择和配准参数与官网输出一致；没有取得官网原始 JSON 的文件哈希，不能再推断两个 JSON 文件逐字节一致。

## 顶层字段

| 字段 | 此文件的值/结构 | 含义 |
|---|---|---|
| `serialNumber` | `H1-8CPZTMY` | 芯片序列号，与运行参数 `--slide` 对应。 |
| `area` | `A1` | 本次使用的 capture area，与 `--area` 对应。 |
| `metadata.version` | `1.0.0` | 对齐文件格式版本，不是 Space Ranger 的版本。 |
| `metadata.slide_design` | `visium_hd_rc1` | Visium HD 芯片设计标识。 |
| `spot_metadata.bin_level` | `8` | 组织选择按 8 µm bins 编码。 |
| `spot_metadata.bitfield_ordering` | `row-major` | 网格按行优先排列。 |
| `spot_metadata.bitfield_encoding` | `base64` | 组织标记是 Base64 编码的二进制位，而非基因序列。 |
| `spot_metadata.bitfield_bit_convention` | `top-left` | 位图的网格方向约定；解码顺序已用官网空间位置表实际验证。 |
| `spot_count` | `702244` | 838 × 838 的完整 8 µm 网格，不是 tissue-bin 数，也不是细胞数。 |
| `oligo` | 长度 117044 的 Base64 文本 | 组织选择位图；1 表示选作组织，0 表示未选作组织。 |
| `transform` | 3 × 3 矩阵 | 与芯片布局一起用于阵列坐标到图像像素坐标的变换。 |
| `slide_layout_file` | Base64 包装的 JSON | 内嵌该 slide 的设计及区域标定。 |
| `checksum` | `b6eaff0a5aa12551f9783abfeb7ef66c` | 对齐输入图像的校验标识；不是上述整个 JSON 的 MD5，也不能当成两张 TIFF 各自的下载 MD5。 |
| `removeImagePages` | `[]` | 没有指定要移除的图像页。 |

## 组织选择位图：oligo

Base64 解码得到 87781 字节，即 702248 位。
其中前 702244 位为网格标记，最后 4 位是 0 填充。

- 完整网格：702244 个 8 µm bins。
- 组织内：393543 个，约 56.0408%。
- 非组织：308701 个。

这个字段名虽然叫 `oligo`，在此 HD 文件内它表达的是组织/非组织掩膜，不是原始 RNA 或探针序列。

与既有输出的标签核验：

| 对比对象 | 分辨率 | 网格位置总数 | 输出组织位置数 | 与此 JSON 的标签差异数 |
|---|---:|---:|---:|---:|
| 官网 3.0.0 | 8 µm | 702244 | 393543 | 0 |
| 官网 3.0.0 | 2 µm | 11222500 | 6296688 | 0 |
| 用户自动运行 4.0.1 | 8 µm | 702244 | 397587 | 5502 |
| 用户自动运行 4.0.1 | 2 µm | 11222500 | 6336843 | 80663 |

对于 2 µm 检查，使用 8 µm mask 中 `[array_row // 4, array_col // 4]` 对应的标签。
这项核验说明官网的组织选择与此文件一致；不意味着改用此文件后所有版本的表达矩阵都会完全相同。

## transform 与 slide_layout_file

顶层 `transform` 为：

```text
[[ 0.2139317019361344,     -0.00023583183650092922, 351.5297540823815  ],
 [ 0.00041276239096999434,  0.2138904425868365,     761.3386633021917  ],
 [ 0.00000006313211048862358, 0.00000003841764791746863, 1.0           ]]
```

它是空间配准变换，不是图像本身；会与区域标定、spot 的实际间距/中心位置组合。
矩阵含缩放、旋转/倾斜、平移及微小的透视项，不能直接把其中的平移值当作显微镜图像位移。

内嵌 `slide_layout_file` 解码后包含：

- `slide_uid = H1-8CPZTMY`。
- `file_format = 1.0`，`aligner_version = 1.0`。
- `slide_design = visium_hd_rc1`。
- `input_hash = f90c4491f4d5d7b6ad662bc9cdf4d205`，属于布局元数据，不能当作 FASTQ 或图像 MD5。
- `capture_areas` 同时列出 A1、D1 标定；本文件通过 `area = A1` 选择 A1，出现 D1 不表示本次分析用了两块组织。

将 A1 的 9 个标定数字按行恢复为矩阵 L，顶层矩阵记作 T，2 µm spot 中心网格转换记作 C：

```text
C = [[2, 0, 1],
     [0, 2, 1],
     [0, 0, 1]]
```

计算 T × L × C 并归一化后，得到官网 feature_slice.h5 中的 `spot_colrow_to_cytassist_colrow`：

```text
[[ 0.4278609186109019,      -0.000026543962897698465, 360.9474240416034 ],
 [ 0.00038127752892269564,   0.42776736133374704,     770.0102013225131 ],
 [ 0.00000013014881468549344, 0.0000000756904126615235, 1.0             ]]
```

## 它没有保存什么

- 没有基因表达计数、聚类结果或细胞类型。
- 没有两张原始图像本身。
- 没有 `cytAssistInfo` 字段，因此没有以该字段保存的 CytAssist ↔ H&E 显微镜图像手工配准矩阵。
- 没有两张 TIFF 各自的原始文件下载 MD5。

官方公开的 `LoupeParser` 将 HD 组织位图、spot 配准和 `cytAssistInfo` 内的图像配准分开处理。
因此此前笼统的“官网做了手工配准”应细化为：官方输出使用了此类文件的阵列定位/组织选择；不能据此声称 H&E ↔ CytAssist 也做了手工配准。
同样，读取此文件不能断定所有标记最初都是由人手工修改，而非自动结果后直接导出。

## checksum 与普通 MD5 的区别

官方 Python 代码通过 `tiffer checksum <image>` 计算对齐输入图像校验值。
不能仅因为 `checksum` 是 32 位十六进制，就把它当成整个 JSON 文件或两张 TIFF 各自的 `md5sum`。
本次只读取用户授权的 JSON，没有为验证这个字段额外读取 `/data` 下的其他文件或执行 `/data` 中的软件。

## 实际使用

参数可以直接指定原文件，无需因为官网命令中的 basename 不同而改名：

```bash
--loupe-alignment=/data/04_ST_introduction/data/spaceranger/Visium_HD_Mouse_Brain_alignment_file.json
```

未执行重跑。即使增加该参数，4.0.1 与官网 3.0.0 的 barcode/UMI 算法差异仍可能存在。

## 官方字段解释依据

- https://raw.githubusercontent.com/10XGenomics/spaceranger/main/lib/python/cellranger/spatial/loupe_util.py
- https://raw.githubusercontent.com/10XGenomics/spaceranger/main/lib/python/cellranger/spatial/tiffer.py

公开 main 分支代码用于解释字段；不将其当作已安装 3.0.0/4.0.1 的逐字节代码快照。
组织掩膜和变换的匹配结论独立来自本地实际文件数值比较。
