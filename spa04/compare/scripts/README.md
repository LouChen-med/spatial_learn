# scripts

- `compare_spaceranger_results.py`：从官方/本地 HTML 读取命令、从 HDF5 读取版本，比较 QC、2/8/16 µm raw/filtered 矩阵、空间 barcode、基因 pseudobulk 以及 tissue mask。
- 从项目根目录运行：`/home/lcpyspa/miniconda3/bin/conda run --no-capture-output -n cosmx-spatial python spa04/compare/scripts/compare_spaceranger_results.py`。
- 默认结果位置为 `../results/`；`--results-dir` 可指定 `spa04/compare/` 内的结果目录。
- 使用 sparse HDF5 数据分块累计，未将巨大的 gene×bin 矩阵转成稠密矩阵；图像按 array row/column 绘制，不是显微镜图像方向。
