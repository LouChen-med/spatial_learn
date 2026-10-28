# Visium HD Mouse Brain：本地与官网结果比较

## 文件结构

- `AGENTS.md`：任务背景、环境、比较边界与进度。
- `tutorial.md`：针对本数据集的 Space Ranger 命令说明及严格复现条件。
- `scripts/`：只读比较脚本，详见目录内 README。
- `results/`：比较报告、CSV 和图片，详见目录内 README。

## 使用

先看 `results/assessment.md` 和 `results/tissue_mask_comparison.png`。
复算命令如下；使用现有环境，不安装软件、不重跑 Space Ranger：

```bash
cd /home/lcpyspa/spa
/home/lcpyspa/miniconda3/bin/conda run --no-capture-output -n cosmx-spatial \
  python spa04/compare/scripts/compare_spaceranger_results.py
```

输入只读取 `spa04/Visium_HD_Mouse_Brain/outs` 和 `spa04/out_10x`。
依赖为 `numpy`、`h5py`、`pyarrow`、`matplotlib`；已经在上述环境内验证可导入。
没有修改 Conda 环境或 Git，因而没有新增环境/Git 修改命令记录。
