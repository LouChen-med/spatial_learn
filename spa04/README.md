# spa04

- `ST_introduction_lc.ipynb`：空间组学入门学习 Notebook，保留现有代码与输出。
- `ST_introduction_lc.py`：Notebook 的逐单元导出及语法修正版，以 `#%%` 分隔单元；Shell 命令通过标准库 `subprocess` 调用，示例日志和说明文字保留为注释。可由普通 Python 解析，也适合在 VS Code 的 Jupyter 交互窗口中按单元运行；未修改原始 Notebook。
- `.ipynb_checkpoints/`：本地 Jupyter 自动检查点，不上传 Git。

在 Jupyter 或 VS Code 中打开 Notebook，检查依赖与数据路径后按单元格顺序运行。本次未重新运行。

## Python 脚本使用

- 先运行包含 `import subprocess` 的首个代码单元，再按需运行后续单元。
- 每个 Shell 单元使用独立 Bash 进程；`cd`、`export` 和 `source` 的效果不会传递给其他单元或 Python 主进程。
- 原教程中的 `/data/...`、`/saw/...`、`/home/chenao/...` 和 `/path/to/...` 仍是示例路径，需按实际数据位置调整；当前项目操作应限制在 `/home/lcpyspa/spa` 内。
- 安装示例保留了修改 `~/.bashrc` 的命令，仅用于展示原教程；按项目约束，不应直接执行这些写入，应改为当前进程的环境配置。
- 不要直接整文件运行：其中包含软件下载、大数据下载和高资源分析流程，且需要 Scanpy、Stereopy 及各平台命令行软件。
- 2026-10-27：修复 Python/Shell 混写、裸日志及说明文字、Shell 续行符后的空格；仅做静态语法验证，不运行下载或分析，不修改 Conda 环境或 Git。

## Visium HD 输入与输出

- `ST_introduction_lc.py` 的 `spaceranger count` 单元使用 `DATA_DIR` 统一指定输入数据目录。目前沿用 `/data/04_ST_introduction/data/spaceranger`，未访问或确认该项目外目录；项目内暂未发现参考目录、FASTQ、探针集和图像输入。
- `--transcriptome` 指向解压后的完整 10x 参考目录，`--fastqs` 指向解压后的 FASTQ 目录；另外三个文件参数分别指定探针集 CSV、CytAssist 图像和显微镜图像。图片文件名沿用原教程，运行前需核对实际来源与图像类型。
- `OUTPUT_DIR=/home/lcpyspa/spa/spa04/Visium_HD_Mouse_Brain`，通过 `--output-dir` 显式指定；`cd` 使用其父目录 `spa04/`，不会额外嵌套一层同名目录。
- `Visium_HD_Mouse_Brain/`：Space Ranger 运行及结果目录；本次检查为空。完成后主要输出在 `outs/`，例如 `outs/web_summary.html` 和 `outs/metrics_summary.csv`。不要把原始输入混放在这个运行目录内。
- 结果查看单元已同步到该输出路径。分析命令移除了后台 `&`，让 Python 等待任务完成并报告返回码。
- 2026-10-27：仅修改命令、说明及结果查看路径，未启动分析、下载输入、调整 Conda 环境或执行 Git 操作。

参数与输出格式依据 10x 官方文档：

- https://www.10xgenomics.com/support/software/space-ranger/4.0/analysis/running-pipelines/command-line-arguments
- https://www.10xgenomics.com/support/software/space-ranger/4.0/analysis/outputs/output-overview
