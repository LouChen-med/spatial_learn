# %% [markdown]
# # 文献复现
#
# **High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer**
#
# **结直肠癌免疫细胞群的高清空间转录组学分析**

# %%
from datetime import date
print("Author: Biomamba_pudding")
print("Date:", date.today())

# %% [markdown]
# # 零、欢迎关注**Biomamba 生信基地**
#
# ## Biomamba，理学博士，此公众号旨在分享生信知识及科研经验与体会，目前在B站的同名账号中发布的手把手教你做单细胞测序数据分析系列视频点击量已经超过45.8W。同名微信公众号拥有超过五万的关注用户以及超过五百万的访问量。如果需要数据分析指导、生信热点全文复现、自测/数据库数据个性化分析辅导、常态化实验学习，欢迎联系客服微信[Biomamba_yunying]
#
# ![./img/Rmarkdown开头2025.jpg](http://biomamba.com:34080/free.download/file.back/biomamba/Rmarkdown%E5%BC%80%E5%A4%B42025.jpg)

# %% [markdown]
# # 一、写在前面
#
# ## 本教程基于Linux中的jupyter环境演示，计算资源不足的同学可参考：  
# [生信分析为什么要使用服务器？](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491007&idx=1&sn=8619c2bad23791b3252a1c4cb1831931&chksm=9b3f64efac48edf9c02f84959aeee9b12480736fce70c5ddb8b41d79cfc22ce590ea2cb60fba&token=1760636281&lang=zh_CN#rd)\
# [足够支持你完成硕博生涯的生信环境](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511780&idx=1&sn=3acc0874f7d2fa75f1abbe7a9a9b4069&chksm=9b3cd7b4ac4b5ea27d4d53964e4a4a9fa8c806a19103729dd4ff0fa72e621411f46a251d269f&scene=21#wechat_redirect)\
# [配置一个心仪的工作站(硬件+环境配置)](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247519707&idx=1&sn=73dcab26f9b2ed75a954327317d07112&chksm=9b3cf48bac4b7d9df5b5a432fbe70d4f9c13bc8d739aa401732614f022223d1912642a7e4bea&scene=21#wechat_redirect)\
# [独享服务器，生信分析不求人](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503618&idx=1&sn=ca47ca985aa8560f6b38f15f210519e2&chksm=9b3cb652ac4b3f44292560071b7e98e2f6fa6ec03581cd882a4d59b9d7fa063a5b9ac1b923f3&scene=21#wechat_redirect)\
# [为实验室准备一份生物信息学不动产](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247492383&idx=1&sn=b8d912877cb38e8d8c9084d1d9370ded&chksm=9b3c9a4fac4b135926091a8af972c45a17a2fa1f94c704779c60df5f733847036c3797894c1c&scene=21#wechat_redirect)
#
#
# ## Linux和Conda是我们本次课程环境搭建的基础，这两部分我们在前期推送中也有相关内容介绍。如果还没有接触过的同学可以跟随我们[十小时学会Linux](https://mp.weixin.qq.com/s/77IMXW_6QFLpbqSeMLfNJg)、[生信软件管家——conda的安装、使用、卸载](https://mp.weixin.qq.com/s/VwRHks00awHXY1tfoPNX9A)和[利用conda打包/复刻生信环境](https://mp.weixin.qq.com/s/y_NPTabHSDGPw2D-_Izs1A)教程内容提前学习一下。
#
# ## jupyter notebook作为Python的集成环境，可以逐行运行Python代码并返回结果，并且支持markerdown语法编译，可以说是一个学习Python的笔记本。jupyter notebook使用可以参考[码Python神器：jupyter notebook](https://mp.weixin.qq.com/s/bK_Ri16z8YG-1ZyRqkhM9A)。
#
#
# ##  Python作为本次课程**核心编程语言**。它语法简洁、库资源丰富，计算高效，是数据科学、自动化开发等领域的首选，兼顾新手入门与资深者高效需求。因此，同学们在学习这个教程前，需要有一定得Python知识，可以参考[生信Python速查手册](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499921&idx=1&sn=dcd6f95596e3de9bf378d850c217fe78&chksm=9b3cb9c1ac4b30d7bbd44057b181f87541f5904d24f903494ae557e768ae00c208c1b68b45c7&payreadticket=HFPtvn0Hc_3itk_RHA90ykGyhDmPMGdHlmAWpCN5_JoUcbpyL6P82ctGZm0c2MiiBCZcLU0#rd)进行快速学习

# %% [markdown]
# # 二、文献解读

# %% [markdown]
# <div style="border: 2px dotted #000; 
#             padding: 10px;
#             background-color: #eef3f5ff; 
#             border-radius: 5px;
#             height: auto;
#             font-size: 24px; 
#             font-weight: bold; 
#             color: #222; 
#             text-shadow: 0px 0px 1px #fff;">
#
# 题目：<p style="font-size: 40px;"><b>High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer</b></p>
# <p style="font-size: 40px;"><b>结直肠癌免疫细胞群的高清空间转录组学分析</b></p>
#
# 影响因子：29<br>
# 期刊：Nature Genetics<br>
# Doi: [10.1038/s41588-025-02193-3](https://doi.org/10.1038/s41588-025-02193-3).<br>
# 发表时间：2025.6.5<br>
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx1.png)
#
# </div>

# %% [markdown]
# ## 1. 🔎 背景
#
# - 结直肠癌（Colorectal cancer, CRC）是全球发病率 & 死亡率都很高的癌症。虽然已有很多基因组学/单细胞 RNA-seq 的研究，但肿瘤微环境（tumor microenvironment, TME）中免疫细胞在空间上的分布、亚型多样性，以及它们如何与肿瘤细胞／边界区域交互，这些空间上的细节还不够清楚。
#
# - FFPE 样本（石蜡固定石蜡包埋）是临床常见样本类型，但很多空间转录组技术对 FFPE 的兼容性和分辨率有局限。
#
# &emsp;&emsp;全面了解结直肠癌 (CRC) 中的细胞行为和对肿瘤微环境 (TME) 的反应仍然难以捉摸。在这里，我们介绍了高清Visium空间转录组学技术 (Visium HD)，并研究了福尔马林固定石蜡包埋的人CRC样本 (n =5)。我们证明了Visium HD的高灵敏度、单细胞尺度分辨率和空间精度，生成了CRC样本的高度精细的全转录组空间轮廓。我们通过与肿瘤和T细胞的相互作用，在不同的空间生态位中发现了具有潜在促肿瘤和抗肿瘤功能的转录组不同的巨噬细胞亚群。 原位基因表达分析验证了我们的发现，并将克隆性扩增的T细胞群体定位在具有抗肿瘤特征的巨噬细胞附近。我们的研究展示了高分辨率空间技术在理解TME中细胞相互作用的能力，并为更大规模的研究铺平了道路，这些研究将揭示CRC生物学的机制和生物标志物，改善诊断和疾病管理策略。
#
# ---

# %% [markdown]
# ## 2. 🛠 方法
#
# ### 2.1 数据与技术
# - **样本**：5 位 CRC 患者的 FFPE 样本 + 邻近正常组织 (NAT)。  
# - **技术**：Visium HD （高分辨率，2–8 µm bins），对比 Visium v2；Xenium 空间原位杂交
# - **参考数据**：单细胞 RNA-seq atlas，用于空间 deconvolution 和 cell type 注释 
#
# ### 2.2 分析流程
# 1. 🔬 **空间测序**：Visium HD 捕获转录本分布。  
# 2. 🎯 **精度验证**：利用 goblet cell marker (CLCA1, FCGBP, MUC2)，定位精度高达 **98–99%**
# 3. 🧩 **Deconvolution**：结合单细胞参考，将空间 bin 注释为具体细胞类型
# 4. 🏰 **边界区分析**：定义肿瘤边缘 ≤50 µm 的区域，重点研究免疫细胞组成与互作
# 5. ✅ **验证**：in situ (Xenium) + 成像 + 核分割校正
#
# ---

# %% [markdown]
# ## 3. 🌟主要发现
#
# ### 3.1 Visium HD 的高精度性能
#
# #### （1）🧩 技术改进
# - **Visium HD 捕获阵列**  
#   - 2 × 2 μm 的条形码方格，连续排列，无空隙  
#   - 提供比 Visium v2 高几个数量级的条形码复杂度与空间分辨率  
# - **输出数据**  
#   - Space Ranger v3.0 输出：
#     - 原始 **2 μm 数据**
#     - **8 μm** 与 **16 μm binning 数据**  
#   - 本文分析主要使用 **8 μm binning 数据**
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx2.png)

# %% [markdown]
# #### （2）📊 性能评估
#
# | 指标 | Visium v2 | Visium HD | 改进效果 |
# |------|-----------|-----------|-----------|
# | 捕获区域 | 55 μm 圆形点，100 μm 间距 | 2 × 2 μm 方格，连续排列 | 分辨率提升数个数量级 |
# | 分辨率 | 聚类较粗糙 | 聚类更细致，能更好对应形态学特征 | **更高分辨率** |
# | UMI 灵敏度 | 高 | 高 | **灵敏度保持一致** |
# | UMI 相关性 (R²) | ~0.82–0.90 | 与 v2 高度相关 | **无灵敏度损失** |
# | gDNA 噪音 | 平均 4.13% | 平均 0.70% | **显著降低** |
# | 过滤噪音后相关性 (R²) | - | ~0.92–0.96 | **数据更可靠** |

# %% [markdown]
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx3.png)
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx4.png)

# %% [markdown]
#
# - **验证实验**  
#   - 标记基因：杯状细胞 marker **CLCA1, FCGBP, MUC2**  
#   - 结果：  
#     - 98–99% 转录本定位在正确区域  
#     - 仅 ~1% 出现在相邻区域  
#   - → 证明 **Visium HD 可达到近似单细胞的空间精确性**
#
# #### （3）✅ 总结
# Visium HD 相较 Visium v2：
# - **分辨率显著提高**
# - **噪音显著减少**
# - **灵敏度保持稳定**
# - **空间精度接近单细胞水平**
#
# ---

# %% [markdown]
# ### 3.2 Visium HD 在单细胞水平上绘制CRC的全景图
#
# #### （1）📊 数据来源
# - **样本**：3 个 CRC 组织切片，使用 **Visium HD** 测序  
# - **单细胞参考图谱**：260,506 个细胞  
#   - 来源：正常与 CRC FFPE 组织 (n = 8)  
#   - 用途：提供细胞类型注释的参考 
#
# #### （2）🔎 初步分析
# - 合并 3 个样本 → 得到 **23 个 cluster**  
# - 聚类细胞类型分为 **9 大类 + 1 未知类群**：  
#   - 🧩 Tumor  
#   - 🧩 Intestinal epithelial  
#   - 🧩 Endothelial  
#   - 🧩 Smooth muscle  
#   - 🧩 T cells  
#   - 🧩 Fibroblasts  
#   - 🧩 B cells  
#   - 🧩 Myeloid  
#   - 🧩 Neuronal  
#   - ❓ Unknown cluster（基因特征未能对应已知类型）  
#
# 👉 结果与组织形态学特征高度一致 (Fig. 3a)   
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx5.png)

# %% [markdown]
# #### （3）🧭 单细胞参考图谱
# - 对 260,506 个细胞做 **差异表达分析 (DGE)**  
# - 使用 **已知标记基因** 进行注释 (Extended Data Fig. 3)  
# - 将其作为 **参考基准**，用于对 Visium HD 数据做 **去卷积 (deconvolution)**  ——spacexr
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx7.png)
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx6.png)

# %% [markdown]
# #### （4）🧩 去卷积结果
# - 得到 **高分辨率的细胞类型空间分布图**  
# - 与预期组织形态相符 (Fig. 3b)  
# - 每个患者样本都有 **主要 & 独特的肿瘤细胞类型**，定位于肿瘤区域 (Extended Data Figs. 4–5)  
#
# ✅ **空间验证 (marker 基因表达一致性)：**  
# - PIGR → goblet cells & enterocytes  
# - CEACAM6 → tumor  
# - COL1A1 → fibroblasts (Fig. 3c) 
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx8.png)
#
# ---

# %% [markdown]
# ### 3.3 🧩 肿瘤边界分析揭示的巨噬细胞亚群
#
# #### （1）🎯 分析思路
# - 研究对象：CRC 组织中 **肿瘤边界 (≤50 μm)** 区域的免疫细胞  
# - 方法：  
#   1. 通过 **spot deconvolution** 标注肿瘤区域  
#   2. 提取所有距离肿瘤 ≤50 μm 的 barcoded 8-μm bins (Fig. 4a)  
#   3. 对比该区域与远离肿瘤区域的 **细胞组成差异**  
#
# #### （2）📊 主要发现
# - **CAF (癌相关成纤维细胞)**：肿瘤边界最丰富的细胞类型  
# - **巨噬细胞**：边界区域最丰富的免疫细胞 (Fig. 4b)  
#   - 验证：C1QC (巨噬细胞标记)，COL1A1 (CAF 标记) (Fig. 4a)  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx9.png)

# %% [markdown]
# #### （3）🧬 巨噬细胞亚群
# 在肿瘤微环境 (TME) 中发现 **两个主要亚群** (Fig. 4c,4d)：  
# - **SELENOP⁺ 巨噬细胞**  
#   - 富集通路：炎症反应、TNF–NF-κB 信号、凋亡、UV 应答  
# - **SPP1⁺ 巨噬细胞**  
#   - 富集通路：凝血、胆固醇稳态、KRAS 信号  
# - 共同特征：EMT 通路和补体系统均上调 (Fig. 5a)  
#
# 👉 两个亚群在空间上主要分布在 **不同的微环境“生态位”**
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx10.png)

# %% [markdown]
# #### （4）🏗️ 巨噬细胞与肿瘤细胞
# - **SPP1⁺ 巨噬细胞邻近的肿瘤细胞** → 上调 **TGFBI, PERP** (Fig. 5b)  
# - **SELENOP⁺ 巨噬细胞邻近的肿瘤细胞** → 上调 **REG1A, REG1B (P1CRC), LCN2 (P2CRC & P5CRC)**  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx11.png)
#
# 📌 **关键发现**：  
# - P1CRC 和 P2CRC 中，REG1A/REG1B 表达的 **杯状细胞 (Goblet cells)** 位于 **肿瘤-正常组织边界**  
# - 这些 goblet cells 可能正在向肿瘤状态转化 (Fig. 5c, 5d)  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx12.png)

# %% [markdown]
# #### （5）🔬 巨噬细胞多样性扩展分析
# - 在全组织范围识别到 **4 个巨噬细胞亚群**：  
#   - **SPP1⁺**：仅限 1 个亚群  
#   - **SELENOP⁺**：分布于 3 个亚群，其中一个低水平表达 **FOLR2**  
# - **FOLR2⁺ 巨噬细胞**：几乎全部 (>98%) 位于 **>50 μm 远离肿瘤** 的区域 (Extended Data Fig. 7)  
#
# ---
#
# #### （6）🔗 细胞通讯 (LIANA 分析, Fig. 5e,f)
# - **SPP1⁺ 巨噬细胞**  
#   - 与肿瘤/T 细胞的主要互作 → **CD44 受体**  
# - **SELENOP⁺ 巨噬细胞**  
#   - 与肿瘤/T 细胞的互作 → 主要影响 **促肿瘤代谢通路 & 免疫刺激**  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx13.png)
#
# ---

# %% [markdown]
# ### 3.4 🧬 CRC TME 中的 T 细胞克隆扩增生态位
#
# #### （1）🎯 研究背景
# - **T 细胞浸润**与 CRC 预后密切相关  
# - **功能状态 & 克隆扩增** 对抗肿瘤作用至关重要  
#
# #### （2）🔍 方法与分析流程
# 1. **肿瘤边界 T 细胞检测**
#    - 使用 **8 × 8 μm bins** → 许多 T 细胞被标记为 **doublet / rejected**，难以定位 (Fig. 6a)  
#    - 采用 **2 × 2 μm bins + 核分割 (segmented nuclear polygons)** → 构建 **gene-by-nucleus UMI 矩阵** → 精准定位 (Fig. 6b)  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx14.png)
#
# 2. **细胞类型标记验证**
#    - T 细胞标记：**TRAC, CD3E**  
#    - 其他细胞标记：PECAM1 (内皮)，IGKC (浆细胞)，COL1A1 (CAF)，SPP1 / SELENOP (巨噬细胞)，CEACAM5 (肿瘤)  
#    - 定位到 **肿瘤边界的 CD4⁺ & CD8⁺ T 细胞** (Fig. 6c)  
#    - 正常组织中未发现 → 提示其 **潜在抗肿瘤作用**  
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx15.png)
#
# ---

# %% [markdown]
# ## 3.5 📝 核心结论
#
# - **Visium HD** 提供单细胞级别高保真空间转录组  
# - **肿瘤边界 T 细胞克隆生态位**
#   - 与巨噬细胞、CAF、肿瘤细胞形成复杂互作  
#   - 可能参与抗肿瘤免疫应答或免疫调控  
# - **免疫细胞空间异质性**
#   - 巨噬细胞亚群 SELENOP⁺ / SPP1⁺ 展现不同功能与局部化  
#   - CAF 与肿瘤细胞协同塑造微环境  
# - **技术互补性**
#   - Visium HD → 全转录组覆盖  
#   - Xenium → 高灵敏度、精细空间定位  
#
# ---

# %% [markdown]
# # 三、代码复现

# %% [markdown]
# ## Fig1：使用Visium HD对CRC和NAT样本进行分析
#
# CRC：colorectal cancer 结直肠癌
#
# NAT: normal adjacent tissue 正常癌旁组织
#
# ![](http://biomamba.com:34080/free.download/file.back/biomamba/chenao/reproduction/wx16.png)

# %% [markdown]
# ### 下载所需数据
#
# 根据文章信息我们可以在 NCBI [GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318) 或者 [10X官方数据库](https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc)下载所需数据，这里我们已经下载好，读取相应数据

# %%
# 导入分析所需的库

import os, re, json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt
import pymde
from PIL import Image

plt.rcdefaults()

# %% [markdown]
# ### 读取单细胞转录组参考数据

# %%
adata_merged_sc = sc.read_10x_h5('/home/chenao/jupyter_home/downstream/process/reproduction/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5')

meta_sc = pd.read_csv('/home/chenao/jupyter_home/downstream/process/reproduction/SingleCell_Metadata_merged.csv')

# %%
adata_merged_sc

# %%
meta_sc

# %%
meta_sc.set_index('Barcode', inplace=True)

# %%
print(meta_sc.index.isin(adata_merged_sc.obs_names).sum(), len(meta_sc)) 

# %%
adata_merged_sc.obs = adata_merged_sc.obs.join(meta_sc)

# %%
adata_merged_sc.obs

# %%
adata_merged_sc.obsm['X_umap'] = adata_merged_sc.obs[['UMAP1', 'UMAP2']].to_numpy()

# %%
sc.pl.umap(adata_merged_sc, color=['Level1'])

# %%
from sklearn.neighbors import KernelDensity
def enhance_umap(adata_merged, color='celltype', figsize=(6,6)):
    fig, ax = plt.subplots(figsize=figsize)

    # 绘制基础 UMAP 点图，不显示 legend、坐标轴框
    sc.pl.umap(
        adata_merged,
        color=color,
        legend_loc=None,
        frameon=False,
        size=15,
        ax=ax,
        add_outline=True,
        show=False
    )
    for cell_type in adata_merged.obs[color].unique():
        # 取出该细胞类型的 UMAP 坐标
        mask = adata_merged.obs[color] == cell_type
        x = adata_merged.obsm["X_umap"][mask, 0]
        y = adata_merged.obsm["X_umap"][mask, 1]

        # 创建网格，用于 KDE
        x_min, x_max = x.min() - 0.5, x.max() + 0.5
        y_min, y_max = y.min() - 0.5, y.max() + 0.5
        xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        xy_train = np.vstack([x, y]).T
        xy_test = np.vstack([xx.ravel(), yy.ravel()]).T

        # 核密度估计 (Kernel Density Estimation)
        kde = KernelDensity(bandwidth=0.5, metric='euclidean')
        kde.fit(xy_train)

        Z = np.exp(kde.score_samples(xy_test))
        Z = Z.reshape(xx.shape)

        # 绘制轮廓线（只绘制密度最高 10% 的边界）
        plt.contour(
            xx, yy, Z,
            levels=[Z.max() * 0.1],
            colors='gray',
            linestyles='--',
            linewidths=1.2,
            alpha=0.75
        )
    cell_type_counts = adata_merged.obs[color].value_counts(normalize=True) * 100
    for cell_type, percentage in cell_type_counts.items():
        mask = adata_merged.obs[color] == cell_type
        x = np.median(adata_merged.obsm["X_umap"][mask, 0])
        y = np.median(adata_merged.obsm["X_umap"][mask, 1])

        ax.text(
            x, y, f"{cell_type}",  # 也可以加百分比 f"{cell_type} ({percentage:.1f}%)"
            fontsize=10,
            ha="center",
            color="black",
            bbox=dict(
                facecolor='white',  # 白色背景
                edgecolor='none',
                alpha=0.65,
                boxstyle="round,pad=0.4"  # 圆角
            )
        )
    ax = plt.gca()
    ax.set_xlabel('UMAP1', fontsize=10)
    ax.set_ylabel('UMAP2', fontsize=10)

    # 获取坐标轴范围
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_range = x_max - x_min
    y_range = y_max - y_min
    aspect_ratio = (fig.get_size_inches()[0] / fig.get_size_inches()[1]) * (y_range / x_range)

    # 箭头长度
    dx = x_range * 0.1 * aspect_ratio
    dy = y_range * 0.1 * aspect_ratio
    head_width = dy * 0.1
    head_length = dx * 0.1

    # X轴箭头
    ax.arrow(x_min, y_min, dx, 0, head_width=head_width, head_length=head_length, fc='black', ec='black')
    # Y轴箭头
    ax.arrow(x_min, y_min, 0, dy, head_width=head_width, head_length=head_length, fc='black', ec='black')

    # 箭头标签
    label_offset = dx * 0.15
    ax.text(x_min + dx/2, y_min - label_offset, 'UMAP1', ha='center', va='top', fontsize=10)
    ax.text(x_min - label_offset, y_min + dy/2, 'UMAP2', rotation=90, ha='right', va='center', fontsize=10)


# %%
enhance_umap(adata_merged_sc, color='Level1', figsize=(6,6))

# %% [markdown]
# ### 读取10x Visium HD 空间转录组数据

# %%
from pathlib import Path
# 文件目录中批量读取数据
adata_merged_dct = {}
meta_dct = {}

data_dir = Path("/home/chenao/jupyter_home/downstream/process/reproduction/HD_data")

for sample_dir in data_dir.iterdir():  # 只遍历一级子目录
    if sample_dir.is_dir():
        _s = sample_dir.name  # 取文件夹名作为样本名
        print(f"Reading sample: {_s}")
        _a = sq.read.visium(str(sample_dir), library_id=_s)
        _a.var_names_make_unique()
        _a.obs_names = [f"{_s}-{i}" for i in _a.obs_names]
        adata_merged_dct[_s] = _a

        # 读取 metadata_merged 并加前缀
        meta_file = sample_dir / "Metadata_merged.parquet"
        if meta_file.exists():
            meta_f = pd.read_parquet(meta_file)
            # 给索引加前缀保证唯一
            meta_f.index = [f"{_s}-{i}" for i in meta_f['barcode']]
            meta_dct[_s] = meta_f
        else:
            print(f"No Metadata_merged.parquet found for sample {_s}")

# %%
meta_dct

# %%
adata_merged_merged = sc.concat(adata_merged_dct,label="library_id",uns_merge="unique")

# %%
adata_merged_merged.obs

# %%
meta_merged = pd.concat(meta_dct.values(), axis=0)

# %%
meta_merged

# %%
meta_merged = meta_merged.reindex(adata_merged.obs_names)

# %%
adata_merged.obs = pd.concat([adata_merged.obs, meta_merged], axis=1)

# %%
adata_merged.obs

# %% [markdown]
# ### 读取 Xenium In Situ gene 对象

# %%

# %% [markdown]
# ## Fig3：使用Visium HD对CRC样本进行空间定位显示出高分辨率、准确的转录本定位

# %% [markdown]
# ### 可视化细胞类型

# %%
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id==sample_id],
    color=["DeconvolutionLabel1"],
    size=2,
    library_id=sample_id,
    figsize=(12,8)
)

# %% [markdown]
# 发现输入的img图像是整个组织切片的，但覆盖spot的区域较小，为了方便我们可视化，需要方法覆盖的区域

# %%
#我们可以通过 空间转录组对象 obsm['spatial'] 矩阵获取 X 和 Y 的最大坐标

samples = adata_merged.obs['library_id'].unique()
ranges = {}

for s in samples:
    coords = adata_merged[adata_merged.obs['library_id'] == s].obsm['spatial']
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    
    padding = 50
    x_min = max(0, x_min - padding)
    y_min = max(0, y_min - padding)
    x_max = x_max + padding
    y_max = y_max + padding
    
    # 用 tuple 保存，sq.pl.spatial_scatter 可以直接使用
    ranges[s] = (x_min, y_min, x_max, y_max)

for name, r in ranges.items():
    print(f"{name} X 范围: {r[0]} - {r[2]}, Y 范围: {r[1]} - {r[3]}")

# %%
# 通过控制 library_id 和 crop_coord 两个参数来可视化我们想要展现的图像
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id==sample_id],
    color=["DeconvolutionLabel1"],
    size=2,
    crop_coord=ranges[sample_id],
    library_id=sample_id,
    figsize=(12,8)
)

# %%
sample_id = 'P2CRC'
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id==sample_id],
    color=["DeconvolutionLabel1"],
    size=2,
    crop_coord=ranges[sample_id],
    library_id=sample_id,
    figsize=(12,8)
)

# %%
sample_id = 'P5CRC'
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id==sample_id],
    color=["DeconvolutionLabel1"],
    size=2,
    crop_coord=ranges[sample_id],
    library_id=sample_id,
    figsize=(12,8)
)

# %% [markdown]
# ### 基因表达可视化

# %% [markdown]
# #### 常规标准化流程

# %%
adata_merged_merged.layers['raw_counts']=adata_merged_merged.X.copy()
sc.pp.normalize_total(adata_merged_merged)
sc.pp.log1p(adata_merged_merged)

# %%
import warnings
warnings.filterwarnings("ignore")  # 忽略所有 warning

# %%
gene_plot = ["PIGR","CEACAM6","COL1A1"]
samples = ['P1CRC','P2CRC','P5CRC']

fig, axes = plt.subplots(len(samples), len(gene_plot), figsize=(12, 12), squeeze=False)

for i, sample_id in enumerate(samples):
    adata_merged_sub = adata_merged_merged[adata_merged_merged.obs.library_id == sample_id]
    
    for j, gene in enumerate(gene_plot):
        ax = axes[i, j]
        sq.pl.spatial_scatter(
            adata_merged_sub,
            color=gene,
            cmap='Reds',
            size=3,
            crop_coord=ranges[sample_id],
            library_id=sample_id,
            img_cmap="Greys",
            ax=ax
        )
        if j == 0:
            ax.set_ylabel(sample_id)
        ax.set_title(gene)

plt.tight_layout()
plt.show()

# %%
adata_merged.write("/home/chenao/jupyter_home/downstream/process/reproduction/merge_HD.h5ad")

# %%
adata_merged = sc.read_h5ad('/home/chenao/jupyter_home/downstream/process/reproduction/merge_HD.h5ad')

# %% [markdown]
# ## Fig4：每个CRC区段中肿瘤周边的细胞组成

# %% [markdown]
# ### Fig4a 可视化肿瘤与周边细胞
#
# 作者这里对肿瘤周边组织进行分析，将注释为肿瘤细胞的bins ( 8μm )显示为红色，并定义肿瘤周边50 μ m内的bins显示为蓝色

# %%
# 查看作者对 spot 的类型定义（这里视研究内容所定）
adata_merged_merged.obs['Periphery'].value_counts()

# %%
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

gene_plot = ["Periphery"]
samples = ['P1CRC','P2CRC','P5CRC']

# 定义颜色映射
color_dict = {
    "50 micron": "blue",
    "Tumor": "red",
    "Tissue": "lightgrey"
}

# %%
fig, axes = plt.subplots(len(samples), len(gene_plot), figsize=(12, 12), squeeze=False)

for i, sample_id in enumerate(samples):
    adata_merged_sub = adata_merged_merged[adata_merged_merged.obs.library_id == sample_id].copy()
    
    # 确保 categorical
    adata_merged_sub.obs["Periphery"] = adata_merged_sub.obs["Periphery"].astype("category")
    adata_merged_sub.obs["Periphery"].cat.set_categories(list(color_dict.keys()))
    
    # 构造颜色列表
    colors = [color_dict[cat] for cat in adata_merged_sub.obs["Periphery"].cat.categories]
    cmap = ListedColormap(colors)
    
    for j, gene in enumerate(gene_plot):
        ax = axes[i, j]
        sq.pl.spatial_scatter(
            adata_merged_sub,
            color=gene,
            palette=cmap,
            size=1,
            crop_coord=ranges[sample_id],
            library_id=sample_id,
            img_cmap="gray",
            ax=ax,
            legend_loc=None   # 不要每张图的图例
        )
        if j == 0:
            ax.set_ylabel(sample_id)
        ax.set_title(gene)

# ====== 全局 legend ======
handles = [mpatches.Patch(color=color, label=label) for label, color in color_dict.items()]
fig.legend(handles=handles, loc='upper center', ncol=len(color_dict), frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.95])  # 给 legend 留出空间
plt.show()

# %%
import numpy as np

def get_crop_bins(adata_merged, center_barcode, size_microns=500):
    """
    返回指定中心点和边长的 crop 坐标以及落在区域内的 spot 数量

    Parameters
    ----------
    adata_merged : AnnData
        包含空间坐标的 AnnData 对象，obsm['spatial'] 已经存放坐标
    center_barcode : str
        用作方框中心的 spot barcode
    size_microns : float
        方框边长，单位 μm

    Returns
    -------
    crop_coord : tuple
        (xmin, ymin, xmax, ymax)
    bins_inside : list
        落在方框内的 spot barcodes
    """
    # 获取中心坐标
    if center_barcode not in adata_merged.obs_names:
        raise ValueError(f"{center_barcode} not found in adata_merged.obs_names")
    
    idx = adata_merged.obs_names.get_loc(center_barcode)
    cx, cy = adata_merged.obsm['spatial'][idx]

    half = size_microns
    xmin, xmax = cx - half, cx + half
    ymin, ymax = cy - half, cy + half

    # 找到落在方框内的 spots
    coords = adata_merged.obsm['spatial']
    inside_mask = (coords[:, 0] >= xmin) & (coords[:, 0] <= xmax) & \
                  (coords[:, 1] >= ymin) & (coords[:, 1] <= ymax)
    bins_inside = adata_merged.obs_names[inside_mask].tolist()

    return (xmin, ymin, xmax, ymax), bins_inside

# 示例使用：
center = "P1CRC-s_008um_00556_00532-1"
crop_coord, bins = get_crop_bins(adata_merged[adata_merged.obs.library_id=='P1CRC'], center, size_microns=1000)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")

# %%
adata_merged_merged = adata_merged_merged[~adata_merged_merged.obs['DeconvolutionLabel1'].isna()].copy()

# %%
adata_merged_merged

# %%
import matplotlib.colors as mcolors
adata_merged_merged.obs["Periphery"] = adata_merged_merged.obs["Periphery"].astype("category")
adata_merged_merged.obs["Periphery"].cat.set_categories(list(color_dict.keys()))

gray_red_cmap = mcolors.LinearSegmentedColormap.from_list("gray_red", ["lightgray", "red"])

# 直接画图
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id == 'P1CRC'],
    color=["Periphery","C1QC","COL1A1"],
    library_id="P1CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    palette=cmap,
    shape='square',
    cmap=gray_red_cmap, 
    size=1, 
    img=None
)

# %%
# 样本P2CRC：
center = "P2CRC-s_008um_00587_00568-1"
crop_coord, bins = get_crop_bins(adata_merged_merged[adata_merged_merged.obs.library_id=='P2CRC'], center, size_microns=1000)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")

# 直接画图
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id == 'P2CRC'],
    color=["Periphery","C1QC","COL1A1"],
    library_id="P2CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    palette=cmap,
    shape='square',
    cmap=gray_red_cmap, 
    size=1, 
    img=None
)

# %%
# 样本P5CRC：
center = "P5CRC-s_008um_00348_00165-1"
crop_coord, bins = get_crop_bins(adata_merged_merged[adata_merged_merged.obs.library_id=='P5CRC'], center, size_microns=1000)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")

# 直接画图
sq.pl.spatial_scatter(
    adata_merged_merged[adata_merged_merged.obs.library_id == 'P5CRC'],
    color=["Periphery","C1QC","COL1A1"],
    library_id="P5CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    palette=cmap,
    shape='square',
    cmap=gray_red_cmap, 
    size=1, 
    img=None
)

# %% [markdown]
# ### Fig4b 计算不同细胞的比例

# %%
# 只保留 singlet
adata_merged_singlet = adata_merged_merged[adata_merged_merged.obs['DeconvolutionClass'] == 'singlet'].copy()

# %%
df = adata_merged_singlet.obs[['DeconvolutionLabel1', 'Periphery', 'library_id']].copy()

# %%
count_table = pd.crosstab(
    index=[df['library_id'], df['Periphery']],  # 行索引：sample + periphery
    columns=df['DeconvolutionLabel1']           # 列索引：细胞类型
)

# %%
prop_table = count_table.div(count_table.sum(axis=1), axis=0) * 100

# %%
prop_table 

# %%
# 先保留你关心的细胞类型
cell_types = [
    "CD4 T cell", "Plasma", "Enterocyte", "Proliferating Immune II", 
    "Pericytes", "Myofibroblast", "Proliferating Fibroblast", 
    "Proliferating Macrophages", "Neutrophil", "Goblet", "Macrophage", 
    "Endothelial", "CAF"
]

# %%
df_prop = prop_table.reset_index()  # 把 multiindex 变成普通列

# %%
df_long = df_prop.melt(
    id_vars=['library_id', 'Periphery'],
    value_vars=cell_types,
    var_name='CellType',
    value_name='Proportion'
)

# %%
df_long = df_long[df_long['Periphery'] != 'Tumor']
df_long['Periphery'] = df_long['Periphery'].cat.remove_unused_categories() 
df_long['Periphery'].unique

# %%
# 自定义颜色和形状
palette = {"50 micron": "blue", "Tissue": "gray"}  # Periphery 对应颜色
markers = {"P1CRC": "o", "P2CRC": "s", "P5CRC": "^"}  # library_id 对应形状

# 绘图
plt.figure(figsize=(12, 8))
sns.scatterplot(
    data=df_long,
    x='Proportion',
    y='CellType',
    hue='Periphery',
    style='library_id',
    palette=palette,
    markers=markers,
    s=100
)
plt.xlabel("Proportion (%)")
plt.ylabel("Cell Type")
plt.title("Cell Type Composition by Periphery and Patient")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # 图例放右边
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Fig4c 巨噬细胞的肿瘤异质性

# %%
from pathlib import Path
# 文件目录中批量读取数据
adata_merged_dct = {}
meta_dct = {}

data_dir = Path("/home/chenao/jupyter_home/downstream/process/reproduction/HD_data")

for sample_dir in data_dir.iterdir():  # 只遍历一级子目录
    if sample_dir.is_dir():
        _s = sample_dir.name  # 取文件夹名作为样本名
        print(f"Reading sample: {_s}")
        _a = sq.read.visium(str(sample_dir), library_id=_s)
        _a.var_names_make_unique()
        _a.obs_names = [f"{_s}-{i}" for i in _a.obs_names]
        adata_merged_dct[_s] = _a

        # 读取 metadata_merged 并加前缀
        meta_file = sample_dir / "Metadata_merged.parquet"
        if meta_file.exists():
            meta_f = pd.read_parquet(meta_file)
            # 给索引加前缀保证唯一
            meta_f.index = [f"{_s}-{i}" for i in meta_f['barcode']]
            meta_dct[_s] = meta_f
        else:
            print(f"No Metadata_merged.parquet found for sample {_s}")

# %%
adata_merged_merged = sc.concat(adata_merged_dct,label="library_id",uns_merge="unique")

# %%
meta_merged = pd.concat(meta_dct.values(), axis=0)
meta_merged = meta_merged.reindex(adata_merged_merged.obs_names)
adata_merged_merged.obs = pd.concat([adata_merged_merged.obs, meta_merged], axis=1)

# %%
# 提取出巨噬细胞 Macrophages

# %%
adata_merged_merged.layers['raw_counts'] = adata_merged_merged.X
adata_merged_merged.write('/home/chenao/jupyter_home/downstream/process/reproduction/HD_rawdata.h5ad')

# %%
# 1. 定义三个筛选条件
# 条件 A: DeconvolutionClass 必须是 'singlet'
is_singlet = (adata_merged_merged.obs['DeconvolutionClass'] == 'singlet')

# 条件 B: DeconvolutionLabel1 必须是 'Macrophage'
is_macrophage = (adata_merged_merged.obs['DeconvolutionLabel1'] == 'Macrophage')

# 条件 C: Periphery 必须是 '50 micron' (即肿瘤微环境 TME 区域)
is_periphery_50micron = (adata_merged_merged.obs['Periphery'] == '50 micron')

# 2. 组合所有三个条件 (使用 & 逻辑与)
# 只有同时满足 A、B、C 的点位才会被选中
combined_condition = is_singlet & is_macrophage & is_periphery_50micron

# %%
tme_macro_adata_merged = adata_merged_merged[combined_condition, :].copy()

# %%
print(f"原始数据点位数: {adata_merged_merged.n_obs}")
print(f"筛选后的 TME Macrophage (singlet, 50 micron) 点位数: {tme_macro_adata_merged .n_obs}")

# %%
tme_macro_adata_merged.X = tme_macro_adata_merged.layers['raw_counts']

# %%
sc.pp.normalize_total(tme_macro_adata_merged, target_sum=1e4)
sc.pp.log1p(tme_macro_adata_merged)
sc.experimental.pp.highly_variable_genes(
    tme_macro_adata_merged, 
    flavor="pearson_residuals",
    layer='raw_counts',
    n_top_genes=3000,
    subset=False,
    inplace=True,
)
sc.pp.scale(tme_macro_adata_merged, max_value=10)
sc.tl.pca(tme_macro_adata_merged, svd_solver='arpack')
sc.pp.neighbors(tme_macro_adata_merged, n_neighbors=10, n_pcs=30)
sc.tl.umap(tme_macro_adata_merged)
sc.tl.leiden(tme_macro_adata_merged, resolution=0.2, key_added='leiden')

# %%
sc.tl.leiden(tme_macro_adata_merged, resolution=0.3, key_added='leiden')

# %%
tme_macro_adata_merged.obs['leiden']

# %%
# Macro标记基因列表
target_genes = [
    'SELENOP', 'STAB1', 'IGKC', 'SLC40A1', 'FGL2', 
    'MS4A6A', 'MPEG1', 'SPP1', 'MMP7', 'CHI3L1', 
    'LPL', 'FN1', 'IL1RN', 'FBP1', 'APOC1'
]

sc.pl.dotplot(
    tme_macro_adata_merged, 
    var_names=target_genes, # 要绘制的基因
    groupby='leiden', # 按照哪个分组（即聚类结果）来计算表达
    # 以下参数用于美化图表
    #standard_scale='var',      # 按变量（基因）进行标准化，使不同基因间表达水平可比
    cmap='viridis',
    figsize=(12, 6),
    dot_max=0.8,               # 限制气泡的最大尺寸
    title="Key Marker Expression in TME Macrophage Subsets"
)

plt.show() # 显示图表

# %%
tme_macro_adata_merged.obs['MacrophageSubtype']

# %%
tme_macro_adata_merged.obs['Combined_ID'] = (
    tme_macro_adata_merged.obs['library_id'].astype(str) + 
    ' - ' + 
    tme_macro_adata_merged.obs['MacrophageSubtype'].astype(str)
)

# %%
tme_macro_adata_merged.obs['Combined_ID']

# %%
desired_order = [
    'P1CRC - Macrophage-SELENOP+',
    'P2CRC - Macrophage-SELENOP+',
    'P5CRC - Macrophage-SELENOP+',
    'P1CRC - Macrophage-SPP1+',
    'P2CRC - Macrophage-SPP1+',
]
tme_macro_adata_merged.obs['Combined_ID'] = pd.Categorical(
    tme_macro_adata_merged.obs['Combined_ID'],
    categories=desired_order,
    ordered=True
)

# %%
sc.pl.dotplot(
    tme_macro_adata_merged, 
    var_names=target_genes, 
    groupby='Combined_ID', 
    swap_axes=True, 
    cmap='viridis',
    figsize=(6, 9), 
    dot_max=0.8
)

# %% [markdown]
# ### Fig4d 巨噬细胞核密度图

# %%
# 定义目标亚型和 TME 区域
TARGET_SUBTYPE_1 = 'Macrophage-SELENOP+'
TARGET_SUBTYPE_2 = 'Macrophage-SPP1+'
TME_REGION = '50 micron'

# --- 1. 筛选坐标 ---

# 筛选条件 1: SELENOP+ 且位于 TME 区域
condition_sel = (
    (adata_merged.obs['MacrophageSubtype'] == TARGET_SUBTYPE_1) & 
    (adata_merged.obs['Periphery'] == TME_REGION)
)
sel_coords = adata_merged.obsm['spatial'][condition_sel.values]

# 筛选条件 2: SPP1+ 且位于 TME 区域
condition_spp = (
    (adata_merged.obs['MacrophageSubtype'] == TARGET_SUBTYPE_2) & 
    (adata_merged.obs['Periphery'] == TME_REGION)
)
spp_coords = adata_merged.obsm['spatial'][condition_spp.values]

print(f"SELENOP+ TME 点位数: {len(sel_coords)}")
print(f"SPP1+ TME 点位数: {len(spp_coords)}")


# %%
def plot_macrophage_density_with_sq_bg(adata_subset: sc.AnnData, title_prefix: str):
    """
    对给定的 AnnData 子集绘制 Macrophage 亚型的空间 KDE 密度图，
    背景使用 Squidpy 的 sq.pl.spatial_scatter 绘制。
    """
    
    TARGET_SUBTYPE_1 = 'Macrophage-SELENOP+'
    TARGET_SUBTYPE_2 = 'Macrophage-SPP1+'
    TME_REGION = '50 micron'
    
    # 确保子集中包含 TME 和 DeconvolutionLabel1
    if 'Periphery' not in adata_subset.obs.columns or 'MacrophageSubtype' not in adata_subset.obs.columns:
        print(f"[{title_prefix}] 警告: 缺少关键列。跳过绘图。")
        return

    # --- 1. 筛选坐标 ---
    
    # SELENOP+ 且位于 TME 区域
    condition_sel = (
        (adata_subset.obs['MacrophageSubtype'] == TARGET_SUBTYPE_1) & 
        (adata_subset.obs['Periphery'] == TME_REGION)
    )
    sel_coords = adata_subset.obsm['spatial'][condition_sel.values]

    # SPP1+ 且位于 TME 区域
    condition_spp = (
        (adata_subset.obs['MacrophageSubtype'] == TARGET_SUBTYPE_2) & 
        (adata_subset.obs['Periphery'] == TME_REGION)
    )
    spp_coords = adata_subset.obsm['spatial'][condition_spp.values]

    # --- 2. 绘制 KDE 密度图 ---
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 绘图逻辑
    def plot_single_kde_sq_bg(ax, coords, subtype_name, cmap):
        
        # 绘制背景点位：使用 Squidpy 绘制，以 'Periphery' 或 'library_id' 为颜色（可选）
        # 这里使用 'Periphery' 作为颜色，以便能区分 TME 和非 TME 区域（即使我们只关注 TME 区域的点位）
        sq.pl.spatial_scatter(
            adata_subset,
            color='Periphery', 
            alpha=0.5,
            title="", # 标题留给 seaborn 绘制
            ax=ax,
            library_id=lib_id
        )
        
        if len(coords) > 0:
            # 绘制 KDE 热图 (在相同的 Axes 上)
            sns.kdeplot(
                x=coords[:, 0], 
                y=coords[:, 1], 
                cmap=cmap,    
                fill=True, 
                alpha=0.7, 
                ax=ax,
                bw_adjust=0.5
            )
        
        ax.set_title(f"{title_prefix}: {subtype_name} Density in TME")
        ax.set_xlabel("Spatial X Coordinate")
        ax.set_ylabel("Spatial Y Coordinate")
        ax.set_aspect('equal', adjustable='box')
        ax.invert_yaxis()
    
    # 调用绘制 SELENOP+
    plot_single_kde_sq_bg(ax1, sel_coords, "SELENOP+", "YlOrRd")
    
    # 调用绘制 SPP1+
    plot_single_kde_sq_bg(ax2, spp_coords, "SPP1+", "YlGn")

    plt.tight_layout()
    plt.show()


# %%
adata_merged.obs

# %%
# 假设 adata_merged 是您包含 P1CRC, P2CRC, P3CRC 等数据的聚合 AnnData 对象

target_libraries = ['P1CRC', 'P2CRC'] 

# 检查 'library_id' 列是否存在
if 'library_id' not in adata_merged.obs.columns:
    print("错误: adata_merged.obs 中缺少 'library_id' 列。请检查您的元数据。")
else:
    for lib_id in target_libraries:
        print(f"\n--- 正在处理: {lib_id} ---")
        
        # 筛选出当前 library_id 的 AnnData 子集
        condition_lib = adata_merged.obs['library_id'] == lib_id
        
        # 检查是否有数据
        if condition_lib.sum() == 0:
            print(f"未在 adata_merged 中找到 library_id 为 '{lib_id}' 的点位。")
            continue
            
        # 提取子集
        # 注意：使用 .copy() 以确保操作的独立性
        adata_current = adata_merged[condition_lib, :].copy()
        
        # 调用绘图函数
        plot_macrophage_density(adata_current, title_prefix=lib_id)

    print("\n所有指定的 library_id 绘图完成。")

# %% [markdown]
# ## Fig5: TME中两种巨噬细胞亚群的鉴定与定位

# %% [markdown]
# ### Fig5a 两种巨噬细胞的富集差异

# %%
tme_macro_adata

# %%
# 确保 'MacrophageSubtype' 是分类变量
tme_macro_adata.obs['MacrophageSubtype'] = tme_macro_adata.obs['MacrophageSubtype'].astype('category')

# 1. 计算差异表达基因
# method='wilcoxon' 是单细胞数据中常用的非参数检验方法
# rank_genes_groups 会将结果存储在 adata.uns['rank_genes_groups'] 中
sc.tl.rank_genes_groups(
    tme_macro_adata, 
    groupby='MacrophageSubtype', 
    method='wilcoxon', 
    key_added='mac_subtype_markers'
)

# %%
# 假设 rank_results 和 target_group 已定义
rank_key = 'mac_subtype_markers'
target_group = 'Macrophage-SPP1+' 

# --- 1. 直接访问 AnnData.uns 中的结果 (保持不变) ---
rank_results = tme_macro_adata.uns[rank_key]
groups = rank_results['names'].dtype.names # 获取所有分组的名称

# --- 2. 优化：逐个分组提取数据，避免 NumPy 扁平化带来的长度问题 ---

data_list = []
for group_name in groups:
    # 提取当前组的基因名称和分数
    group_names = rank_results['names'][group_name]
    group_scores = rank_results['scores'][group_name]
    
    # 将结果转换为 DataFrame 格式的字典
    df_group = pd.DataFrame({
        'group': group_name,
        'names': group_names,
        'scores': group_scores
    })
    data_list.append(df_group)

# 合并所有组的结果到一个 DataFrame
results_df = pd.concat(data_list, ignore_index=True)

# %%
results_df

# %%
# --- 3. 筛选目标亚型并准备 GSEA 排序列表 (保持不变) ---

# 筛选出目标组 SPP1+ 的基因及其 Scores
ranked_genes = results_df[results_df['group'] == target_group][['names', 'scores']]
ranked_genes.columns = ['Gene', 'Score']

# 清理可能出现的 NaN/空字符串/字节字符串
# 确保所有基因名称都是可读的字符串
ranked_genes['Gene'] = ranked_genes['Gene'].apply(lambda x: x.decode() if isinstance(x, bytes) else str(x))
ranked_genes = ranked_genes.dropna(subset=['Gene'])
ranked_genes = ranked_genes[ranked_genes['Gene'].str.len() > 0] 
ranked_genes = ranked_genes[ranked_genes['Gene'] != 'nan']

# GSEA 准备：以 Gene 为索引，并按 Score 降序排列
# 注意：如果存在重复基因名，取第一个（即最高分）
ranked_list = ranked_genes.set_index('Gene')['Score'].sort_values(ascending=False)

print(f"✅ GSEA 排序列表准备完成。'{target_group}' 的 Top 5 基因:")
print(ranked_list.head())

# %%
import gseapy as gp
gene_sets_to_use = [
    'MSigDB_Hallmark_2020',  # 核心通路，推荐首先使用
    'GO_Biological_Process_2023'
]

# %%
# 运行 GSEA prerank 分析
# rnk: 排序后的基因列表
# gene_sets: 基因集列表
# processes: 用于并行计算的核数
# permutation_num: 置换测试的次数，越多越精确（默认 1000）
# outdir: 结果输出目录
gsea_results = gp.prerank(
    rnk=ranked_list,             
    gene_sets=gene_sets_to_use,         
    processes=4,                 
    permutation_num=1000,        
    outdir='gsea_output_spp1_macrophage', 
    seed=42,                     # 设置随机种子以保证结果可重复
    format='png'                 # 同时生成富集图
)

print("GSEA 富集分析运行完成。")

# %%
gsea_results.res2d

# %%
# --- 从 rank_genes_groups 结果中安全提取 DataFrame ---
rank_key = 'mac_subtype_markers'
rank_results = tme_macro_adata.uns[rank_key]
groups = rank_results['names'].dtype.names

data_list = []
for group_name_iter in groups:
    df_group = pd.DataFrame({
        'group': group_name_iter,
        'names': rank_results['names'][group_name_iter],
        'scores': rank_results['scores'][group_name_iter]
    })
    data_list.append(df_group)
full_rank_df = pd.concat(data_list, ignore_index=True)

full_rank_df['names'] = full_rank_df['names'].apply(lambda x: x.decode() if isinstance(x, bytes) else str(x))
full_rank_df = full_rank_df.dropna(subset=['names'])
full_rank_df = full_rank_df[full_rank_df['names'].str.len() > 0] 
full_rank_df = full_rank_df[full_rank_df['names'] != 'nan']


# --- 为 SPP1+ 和 SELENOP+ 准备 ranked_list ---
ranked_list_spp1 = full_rank_df[full_rank_df['group'] == 'Macrophage-SPP1+'].set_index('names')['scores'].sort_values(ascending=False)
ranked_list_selenop = full_rank_df[full_rank_df['group'] == 'Macrophage-SELENOP+'].set_index('names')['scores'].sort_values(ascending=False)

# --- 定义基因集 ---
gene_sets_to_use = ['MSigDB_Hallmark_2020'] # 示例图看起来主要来自 Hallmark

# --- 运行 GSEA for SPP1+ ---
print("运行 GSEA for Macrophage-SPP1+")
gsea_res_spp1 = gp.prerank(
    rnk=ranked_list_spp1,             
    gene_sets=gene_sets_to_use,         
    processes=4,                 
    permutation_num=1000,        
    outdir=None, # 不输出文件到磁盘，直接返回结果
    seed=42
).res2d

# --- 运行 GSEA for SELENOP+ ---
print("运行 GSEA for Macrophage-SELENOP+")
gsea_res_selenop = gp.prerank(
    rnk=ranked_list_selenop,             
    gene_sets=gene_sets_to_use,         
    processes=4,                 
    permutation_num=1000,        
    outdir=None, 
    seed=42
).res2d

print("两个亚群的 GSEA 均运行完成。")

# %%
# --- 合并结果 ---
# 选择关键列: Term, NES, FDR q-val
gsea_res_spp1_filtered = gsea_res_spp1[['Term', 'NES', 'FDR q-val']].copy()
gsea_res_spp1_filtered['CellType'] = 'SPP1+'

gsea_res_selenop_filtered = gsea_res_selenop[['Term', 'NES', 'FDR q-val']].copy()
gsea_res_selenop_filtered['CellType'] = 'SELENOP+'

# 将两个 DataFrame 合并
combined_gsea_df = pd.concat([gsea_res_spp1_filtered, gsea_res_selenop_filtered], ignore_index=True)

# --- 筛选显著通路 ---
# 筛选 FDR q-val < 0.25 的通路
significant_combined_gsea = combined_gsea_df[combined_gsea_df['FDR q-val'] < 0.5].copy()

if significant_combined_gsea.empty:
    print("警告: 没有发现显著富集的通路 (FDR q-val < 0.25)。请检查 GSEA 参数或数据。")
else:
    # --- 计算绘图值 ---
    # 按照图示，横坐标是 log(adj. P value)，并且方向由富集亚型决定
    # 我们将使用 -log10(FDR q-val) 乘以 NES 的符号来表示方向和强度

    # 避免 log(0) 的情况，将极小的 p 值替换为更小但非零的值
    # 避免 log(0) 的情况，将极小的 p 值替换为更小但非零的值
    significant_combined_gsea['FDR q-val'] = significant_combined_gsea['FDR q-val'].replace(0, np.finfo(float).eps)

    # 确保 'FDR q-val' 是标准的浮点数类型（float/float64）
    # 这一行是关键修正，用于处理可能的非标准或混合数据类型
    significant_combined_gsea['FDR q-val'] = significant_combined_gsea['FDR q-val'].astype(float) 

    # 现在可以安全地计算 -log10
    significant_combined_gsea['log_adj_p_value'] = -np.log10(significant_combined_gsea['FDR q-val'])
    
    # 将 NES 符号乘以 log_adj_p_value，以区分 SPP1+ (负值) 和 SELENOP+ (正值)
    # 这里的逻辑需要根据您的数据和期望的图示进行调整
    # 示例图中：左侧是 SPP1+ 富集，右侧是 SELENOP+ 富集。
    # 这意味着我们希望 SPP1+ 的值是负的，SELENOP+ 的值是正的。
    # 如果 SPP1+ 的 NES > 0，则它富集在 SPP1+ 侧，对应负值。
    # 如果 SELENOP+ 的 NES > 0，则它富集在 SELENOP+ 侧，对应正值。
    
    significant_combined_gsea['PlotValue'] = np.where(
        (significant_combined_gsea['CellType'] == 'SPP1+'),
        -significant_combined_gsea['log_adj_p_value'] * np.sign(significant_combined_gsea['NES']),
        significant_combined_gsea['log_adj_p_value'] * np.sign(significant_combined_gsea['NES'])
    )

    # 仅保留在某个亚型中显著的通路，并去除重复的通路名称
    # 如果一个通路在两个亚型中都显著，我们只取 NES 绝对值最大的那一个
    significant_combined_gsea['abs_NES'] = np.abs(significant_combined_gsea['NES'])
    
    # 按照图示，每个通路只出现一次
    # 如果一个通路在 SPP1+ 和 SELENOP+ 都显著，需要决定显示哪一个
    # 这里我们简化处理：按 PlotValue 绝对值排序，并取每个 Term 的最高值
    plot_data = significant_combined_gsea.loc[significant_combined_gsea.groupby('Term')['abs_NES'].idxmax()].copy()
    
    # 对 PlotValue 进行排序，以便在图中按顺序显示
    plot_data = plot_data.sort_values(by='PlotValue', ascending=True)

    # 检查通路数量，确保不要太多或太少
    if len(plot_data) > 30: # 限制一下显示的通路数量，避免图过长
        plot_data = pd.concat([plot_data.head(5), plot_data.tail(5)])
    
    # --- 打印用于绘图的数据 ---
    print("\n准备用于绘图的数据 (PlotValue):")
    print(plot_data[['Term', 'CellType', 'NES', 'FDR q-val', 'PlotValue']].head())

# %%
if 'plot_data' in locals() and not plot_data.empty:
    plt.figure(figsize=(10, len(plot_data) * 0.5)) # 根据通路数量调整图的高度
    
    # 定义颜色：SPP1+ 对应绿色，SELENOP+ 对应橙色（或根据 NES 符号）
    # 图中的颜色是基于富集方向的，不是基于 CellType 本身
    # 如果 PlotValue < 0，富集在 SPP1+ 侧 (绿色)
    # 如果 PlotValue > 0，富集在 SELENOP+ 侧 (橙色)
    colors = ['#1f77b4' if val < 0 else '#ff7f0e' for val in plot_data['PlotValue']] # 示例颜色，根据需要调整

    sns.barplot(
        x='PlotValue', 
        y='Term', 
        data=plot_data, 
        palette=colors,
        ax=plt.gca()
    )

    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8) # 零线

    plt.xlabel('log (adj. P value)') # 确保标签与图示一致
    plt.ylabel('') # 移除 Y 轴标签，通路名称即是标签
    plt.title('Differential Pathway Enrichment')

    # 添加 SPP1+ 和 SELENOP+ 箭头和标签
    plt.text(plt.xlim()[0], -1.5, 'SPP1⁺', color='green', ha='left', va='center', fontsize=12, weight='bold') # 调整位置
    plt.text(plt.xlim()[1], -1.5, 'SELENOP⁺', color='orange', ha='right', va='center', fontsize=12, weight='bold') # 调整位置
    plt.arrow(
        plt.xlim()[0] * 0.8, -0.8,  # x, y start
        plt.xlim()[0] * 0.2, 0,  # dx, dy
        head_width=0.3, head_length=0.5, fc='green', ec='green', clip_on=False
    )
    plt.arrow(
        plt.xlim()[1] * 0.2, -0.8,  # x, y start
        plt.xlim()[1] * 0.8, 0,  # dx, dy
        head_width=0.3, head_length=0.5, fc='orange', ec='orange', clip_on=False
    )
    plt.yticks(fontsize=10) # 调整 Y 轴字体大小
    
    plt.tight_layout()
    plt.show()
else:
    print("没有足够的数据来绘制差异通路富集图。")

# %% [markdown]
# ### Fig5b 不同样本巨噬细胞亚型特异基因表达小提琴图

# %%
tme_macro_adata.X = tme_macro_adata.layers['raw_counts']

# %%
# 图 b 中所示的基因列表
genes_to_plot = [
    'LCN2', 'REG1A', 'DMBT1', 'PIGR', 'REG1B', 
    'MME', 'LGR5', 'PERP', 'TGFBI', 'KRT23'
] 

# 检查这些基因是否存在于 AnnData 对象的 .var_names 中
valid_genes = [gene for gene in genes_to_plot if gene in tme_macro_adata.var_names]

if len(valid_genes) != len(genes_to_plot):
    missing_genes = set(genes_to_plot) - set(valid_genes)
    print(f"警告：以下基因在数据中缺失，将跳过绘图: {missing_genes}")

# %%
# 图 b 中所示的基因列表
genes_to_plot = [
    'LCN2', 'REG1A', 'DMBT1', 'PIGR', 'REG1B', 
    'MME', 'LGR5', 'PERP', 'TGFBI', 'KRT23'
] 

# 目标分组变量：现在直接使用 Combined_ID
groupby_col = 'Combined_ID' 

# 检查基因是否存在
valid_genes = [gene for gene in genes_to_plot if gene in tme_macro_adata.var_names]

if not valid_genes:
    print("错误：数据中没有找到任何要绘制的基因。请检查基因名是否正确。")
    # 如果没有有效基因，停止执行
else:
    # --- 2. 绘制小提琴图 ---
    
    # 调整整体画布大小，保持 2x5 布局
    plt.rcParams['figure.figsize'] = (20, 10) 
    fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))
    axes = axes.flatten()

    # 循环绘制每个基因的小提琴图
    for i, gene in enumerate(valid_genes):
        if i < len(axes):
            ax = axes[i]
            
            # 使用 Scanpy 绘制小提琴图，仅按 Combined_ID 分组
            sc.pl.violin(
                tme_macro_adata, 
                keys=gene, 
                groupby=groupby_col, 
                ax=ax,
                stripplot=False, # 移除散点
                show=False,
                rotation=45, # 旋转 X 轴标签，防止 Combined_ID 标签重叠
                # 仅在左侧保留 Y 轴标签
                ylabel='Normalized log expression' if i % 5 == 0 else '' 
            )
            
            ax.set_title(gene, fontsize=14)
            ax.set_xlabel('') # 移除 X 轴标签 (Combined_ID)

            # 移除图例（通常 violin plots 分组过多时不需要图例）
            if ax.get_legend():
                ax.get_legend().remove()

    plt.tight_layout()
    plt.show()

# %% [markdown]
# ### Fig5c Goblet亚型表达气泡图

# %%
# 提取goblet细胞

# 1. 定义三个筛选条件
# 条件 A: DeconvolutionClass 必须是 'singlet'
is_singlet = (adata_merged_merged.obs['DeconvolutionClass'] == 'singlet')

# 条件 B: DeconvolutionLabel1 必须是 'Macrophage'
is_Goblet = (adata_merged_merged.obs['DeconvolutionLabel1'] == 'Goblet')

# 2. 组合所有三个条件 (使用 & 逻辑与)
# 只有同时满足 A、B、C 的点位才会被选中
combined_condition = is_singlet & is_Goblet

# %% [markdown]
# 这里我们直接使用作者提供的分群信息

# %%
adata_merged.obs['GobletSubcluster']

# %%
goblet_adata = adata_merged[adata_merged.obs['GobletSubcluster'].notna()].copy()

# %%
# 图中所示的基因列表
goblet_marker_genes = [
    'REG1B', 'REG1A', 'DMBT1', 'DUOX2', 'OLFM4', 
    'CLCA1', 'TFF3', 'FCGBP', 'MUC2'
]

# 检查基因是否存在于 goblet_adata 中
valid_genes = [gene for gene in goblet_marker_genes if gene in goblet_adata.var_names]

if not valid_genes:
    print("错误：在 goblet_adata 中没有找到任何要绘制的基因。请检查基因名是否正确。")
    # 如果没有有效基因，停止执行
    
# 确保分组列是 categorical 类型（Scanpy 的要求）
groupby_col = 'GobletSubcluster'
if goblet_adata.obs[groupby_col].dtype.name != 'category':
    goblet_adata.obs[groupby_col] = goblet_adata.obs[groupby_col].astype('category')

# %%
sc.pl.dotplot(
    goblet_adata, 
    var_names=valid_genes,  # 基因列表
    groupby=groupby_col,  # 亚群分组（0到6）
    color_map='bwr', # 使用蓝白红渐变色，与图示一致 (Blue-White-Red)
    standard_scale='var', # 对每个基因进行标准化，使得颜色代表相对表达量（Scaled expression）
    figsize=(8, 8),
    swap_axes=True, 
    title='Gene Expression Across Goblet Subclusters',
    # 我们可以通过 dot_max 和 dot_min 来调整点图的大小范围
    # dot_max=0.75, 
    # dot_min=0.01,
    show=False
)

plt.tight_layout()
plt.show()

# %% [markdown]
# ### Fig5d Goblet亚型空间可视化分布

# %%
gene_plot = ["GobletSubcluster"]
samples = ['P1CRC','P2CRC','P5CRC']

# %%
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

goblet_subclusters = [f'Goblet-{i}' for i in range(7)]

# 使用 Matplotlib 的 Set1 调色板（定性颜色，区分度高）
cmap_colors = cm.get_cmap('Set1', len(goblet_subclusters))
colors = [mcolors.rgb2hex(cmap_colors(i)) for i in range(len(goblet_subclusters))]

# 构建颜色字典
color_dict = dict(zip(goblet_subclusters, colors))

print("✅ Goblet 亚群颜色字典 (color_dict):")
print(color_dict)

# %%
fig, axes = plt.subplots(len(gene_plot),len(samples),figsize=(15, 5), squeeze=False)

for i, sample_id in enumerate(samples):
    adata_merged_sub = adata_merged[adata_merged.obs.library_id == sample_id].copy()
    
    # 确保 categorical
    adata_merged_sub.obs["GobletSubcluster"] = adata_merged_sub.obs["GobletSubcluster"].astype("category")
    adata_merged_sub.obs["GobletSubcluster"].cat.set_categories(list(color_dict.keys()))
    
    # 构造颜色列表
    colors = [color_dict[cat] for cat in adata_merged_sub.obs["GobletSubcluster"].cat.categories]
    cmap = ListedColormap(colors)
    
    for j, gene in enumerate(gene_plot):
        ax = axes[j, i]
        sq.pl.spatial_scatter(
            adata_merged_sub,
            color=gene,
            palette=cmap,
            size=1,
            crop_coord=ranges[sample_id],
            library_id=sample_id,
            img_cmap="gray",
            ax=ax,
            legend_loc=None   # 不要每张图的图例
        )
        if j == 0:
            ax.set_ylabel(sample_id)
        ax.set_title(gene)

# ====== 全局 legend ======
handles = [mpatches.Patch(color=color, label=label) for label, color in color_dict.items()]
fig.legend(handles=handles, loc='upper center', ncol=len(color_dict), frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.95])  # 给 legend 留出空间
plt.show()

# %% [markdown]
# ### Fig5e 巨噬细胞与免疫细胞通讯弦图

# %%
adata_merged = sc.read_h5ad('/home/chenao/jupyter_home/downstream/process/reproduction/HD_rawdata.h5ad')
adata_filtered = adata_merged[adata_merged.obs['DeconvolutionLabel1'].notna()].copy()

# %%
adata_filtered.obs

# %%
# 将macro的注释返回整体层次
new_celltype_series = adata_merged.obs['MacrophageSubtype'].astype(str).copy()
decon_label_str = adata_merged.obs['DeconvolutionLabel1'].astype(str)
nan_mask = adata_merged.obs['MacrophageSubtype'].isna()
new_celltype_series[nan_mask] = decon_label_str[nan_mask]
new_celltype_series.replace('nan', pd.NA, inplace=True)
adata_merged.obs['celltype'] = new_celltype_series
adata_merged.obs['celltype'] = adata_merged.obs['celltype'].astype('category')
print(adata_merged.obs['celltype'].value_counts())

# %%
target_labels = [
    'Tumor I', 'Tumor II', 'Tumor III', 'Tumor IV', 'Tumor V',
    'Macrophage-SELENOP+', 'Macrophage-SPP1+', 
    'CD8 T cell', 'CD4 T cell'
]
mask_to_keep = adata_merged.obs['celltype'].isin(target_labels)
adata_merged_filtered = adata_merged[mask_to_keep].copy()

# %%
adata_merged_filtered.obs['celltype'].value_counts()

# %%
# 定义样本 ID 列表
samples = ['P1CRC', 'P2CRC', 'P5CRC']

# 定义细胞类型注释列（上一步创建的整合列）
celltype_col = 'celltype'

# 创建保存目录
output_dir = '/home/chenao/jupyter_home/downstream/process/reproduction/cellphonedb/cellphonedb_input_adata/'
os.makedirs(output_dir, exist_ok=True)
print(f"元数据文件将保存到目录: {output_dir}")

# --- 循环提取、保存和生成元数据文件 ---
for sample_id in samples:
    print(f"\n--- 正在处理样本: {sample_id} ---")

    # 1. 检查样本 ID
    if sample_id not in adata_merged_filtered.obs['library_id'].unique():
        print(f"警告：样本 '{sample_id}' 未在 'library_id' 列中找到，跳过。")
        continue

    # 2. 提取该样本的数据
    adata_sample = adata_merged_filtered[adata_merged_filtered.obs['library_id'] == sample_id].copy()
    num_cells = adata_sample.n_obs
    
    # 3. **保存 AnnData 文件 (.h5ad)**
    output_h5ad_filename = os.path.join(output_dir, f'{sample_id}_filtered.h5ad')
    adata_sample.write_h5ad(output_h5ad_filename)
    print(f"AnnData 文件已保存 ({num_cells} spots): {output_h5ad_filename}")

    # 4. **生成并保存元数据文件 (Meta File)**
    meta_df = pd.DataFrame({
        'Cell': adata_sample.obs_names,  # Spot ID
        'cell_type': adata_sample.obs[celltype_col] # 整合后的细胞类型标签
    })
    
    output_meta_filename = os.path.join(output_dir, f'{sample_id}_meta.txt')
    meta_df.to_csv(output_meta_filename, sep='\t', index=False)
    print(f"Meta 文件已保存: {output_meta_filename}")

# %%
import os
import pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import scanpy as sc
plt.rcdefaults()

# %%
# 1. 定义关键变量
samples = ['P1CRC', 'P2CRC', 'P5CRC']
input_dir = '/home/chenao/jupyter_home/downstream/process/reproduction/cellphonedb/cellphonedb_input_adata' # 您保存 h5ad 和 meta 文件的目录
output_root = '/home/chenao/jupyter_home/downstream/process/reproduction/cellphonedb/cellphonedb_results/'    # 所有结果的主输出目录
cpdb_file_path = '/home/chenao/workspace/download_data/cellphonedb/cellphonedb.zip'        # 默认数据库文件

os.makedirs(output_root, exist_ok=True)
print(f"CellPhoneDB 结果将保存到: {output_root} 下的子目录。")


# --- 2. 循环运行 CellPhoneDB 分析 ---
all_results = {}

for sample_id in samples:
    print(f"\n======== 正在运行样本: {sample_id} 的 CellPhoneDB 分析 ========")
    
    # 构造输入文件路径
    meta_file = os.path.join(input_dir, f'{sample_id}_meta.txt')
    counts_file = os.path.join(input_dir, f'{sample_id}_filtered.h5ad')
    
    # 构造输出路径
    sample_out_path = os.path.join(output_root, sample_id)
    os.makedirs(sample_out_path, exist_ok=True)

    # 检查文件是否存在
    if not os.path.exists(meta_file) or not os.path.exists(counts_file):
        print(f"错误：样本 {sample_id} 的输入文件缺失。请检查路径：{meta_file} 或 {counts_file}")
        continue
    
    try:
        # 调用 CellPhoneDB 分析方法
        cpdb_results = cpdb_statistical_analysis_method.call(
            cpdb_file_path = cpdb_file_path,
            meta_file_path = meta_file,
            counts_file_path = counts_file,
            counts_data = 'hgnc_symbol', # 假设您的基因ID是 HGNC 符号
            output_path = sample_out_path,
            threads = 30
        )
        
        print(f" 样本 {sample_id} 分析完成。结果已保存到 {sample_out_path}")
        # 将结果存储在内存中（可选）
        all_results[sample_id] = cpdb_results 
        
    except Exception as e:
        print(f" 样本 {sample_id} 运行失败。错误信息: {e}")
        # 如果是 CellPhoneDB 数据库路径问题，可能需要手动提供路径

print("\n--- 所有样本 CellPhoneDB 分析运行完毕。---")

# %%
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt

from pathlib import Path
plt.rcdefaults()

# %%
# read in the files
# 1) .h5ad file used for performing CellPhoneDB
DATADIR = Path("/home/chenao/jupyter_home/downstream/process/reproduction/cellphonedb/cellphonedb_results/")

adata = ad.read_h5ad("/home/chenao/jupyter_home/downstream/process/reproduction/cellphonedb/cellphonedb_input_adata/P1CRC_filtered.h5ad")

# 2) output from CellPhoneDB
means = pd.read_csv(DATADIR / "P1CRC" / "statistical_analysis_means_09_26_2025_041655.txt", sep="\t")
pvals = pd.read_csv(DATADIR / "P1CRC" / "statistical_analysis_pvalues_09_26_2025_041655.txt", sep="\t")
decon = pd.read_csv(DATADIR / "P1CRC" / "statistical_analysis_deconvoluted_09_26_2025_041655.txt", sep="\t")

# %%
with pd.option_context('display.max_columns', 100):
    print(pvals.head())

# %%
df = pvals 
melt_id_vars = ['interacting_pair', 'gene_a', 'gene_b']
value_vars_pvals = df.columns[13:].tolist()
pvals_long = df[melt_id_vars + value_vars_pvals].melt(
    id_vars=melt_id_vars,
    var_name='cell_pair',
    value_name='pvalue'
)
print(pvals_long.head())

# %%
df_means = means 

# 1. 确定用于标识交互的元数据列
# 这些列将不会被“熔化”，而是作为重复的标识符保留在新表中
melt_id_vars = ['interacting_pair', 'gene_a', 'gene_b']

# 2. 确定包含均值数据的细胞对列
# 根据我们之前的分析，这些列从索引 13 开始
value_vars_means = df_means.columns[13:].tolist()

# 3. 执行 melt 操作
means_long = df_means[melt_id_vars + value_vars_means].melt(
    id_vars=melt_id_vars,
    var_name='cell_pair',    # 新列名：代表细胞对
    value_name='mean_score' # 新列名：代表均值/交互强度
)

print("Means Long Table Head:")
print(means_long.head())

# %%
plot_data = pd.merge(
    pvals_long, 
    means_long, 
    on=['interacting_pair', 'gene_a', 'gene_b', 'cell_pair'] # 使用所有标识符列进行精确合并
)

# %%
plot_data

# %%
# 您提供的所有目标基因
target_gene_pairs = [
    'APOE', 'SORL1', 'LSR', 'SPP1', 'CD44', 'MMP9', 'MMP7'
]
target_genes = set(target_gene_pairs) 

# 目标细胞对列表
target_cell_pairs = [
    'Macrophage-SPP1+|CD4 T cell', 
    'Macrophage-SPP1+|CD8 T cell', 
    'Macrophage-SPP1+|Tumor II'
]

# --- 细胞对筛选条件 ---
is_target_cellpair = plot_data['cell_pair'].isin(target_cell_pairs)

# --- 基因筛选条件 ---
# 交互作用的任一端 (gene_a 或 gene_b) 包含目标基因
is_gene_a_target = plot_data['gene_a'].isin(target_genes)
is_gene_b_target = plot_data['gene_b'].isin(target_genes)
is_target_gene = is_gene_a_target | is_gene_b_target

# --- 最终的双重筛选 ---
# 必须同时满足细胞对条件 AND 基因条件
filtered_data_final = plot_data[
    is_target_cellpair & is_target_gene
].copy()


if filtered_data_final.empty:
    print("警告：在指定的细胞对和目标基因的交叉筛选下，未找到任何匹配的交互记录。")
    print("请再次核对细胞对和基因名称是否与数据中的命名精确一致。")
else:
    print(f"成功筛选出 {len(filtered_data_final)} 条满足所有条件的交互记录。")
    print("\n最终筛选结果头部（包含目标细胞对和目标基因的交互）：")
    print(filtered_data_final)

# %%
pdata = filtered_data_final.dropna(subset=['gene_b']).copy()

# %%
pdata 

# %%
# 动态调整图形高度以容纳所有交互名称
n_interactions = len(pdata['interacting_pair'].unique())
plt.figure(figsize=(9, max(6, n_interactions * 0.5))) 

sns.scatterplot(
    data=pdata,
    x='cell_pair',             # X 轴：细胞对
    y='interacting_pair',      # Y 轴：交互作用名称
    size='mean_score',         # 点的大小：交互强度
    hue='mean_score',          # 点的颜色：显著性
    palette='magma_r',         # 颜色方案 ('_r'表示反转，让大值（显著）颜色更深)
    sizes=(50, 600),           # 气泡大小范围（根据您的均值范围调整）
    alpha=0.9,
    hue_norm=(0, plot_data['mean_score'].max()) # 颜色范围从 0 开始
)

# 优化图表外观
plt.title('Cell-Cell Interaction Dot Plot (Filtered Data)', fontsize=14)
plt.xlabel('Cell Pair', fontsize=12)
plt.ylabel('Ligand-Receptor Interaction', fontsize=12)

# 旋转 X 轴标签，使其更容易阅读
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)

# 放置图例
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.show()

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 假设 plot_data 已经合并完成，包含 'gene_a', 'cell_pair', 'pvalue', 'mean_score'

# 1. 筛选特定基因 A 的交互作用
# 筛选 'gene_a' 等于 'APOE' 或 'B2M' 的行
genes_to_plot = ['APOE', 'B2M', 'LYZ', 'CD14']
filtered_genes = plot_data[plot_data['gene_a'].isin(genes_to_plot)].copy()

# 2. 筛选特定细胞对
# 确保这里的字符串与您的数据中的 cell_pair 列的值完全一致！
cell_pair_to_plot = 'Macrophage-SELENOP+|CD4 T cell' 
# 注意：如果您的分隔符是下划线或其他符号，请相应调整
filtered_data = filtered_genes[filtered_genes['cell_pair'] == cell_pair_to_plot].copy()

# %%
filtered_data

# %%
kpy.plot_cpdb_heatmap(pvals=pvals, figsize=(5, 5), title="Sum of significant interactions")

# %%
kpy.plot_cpdb(
    adata=adata,
    cell_type1="Macrophage-SELENOP+",
    cell_type2=".",
    means=means,
    pvals=pvals,
    celltype_key="celltype",
    genes=["APOE", "B2M", "LYZ"],
    default_style=False,
    figsize=(6, 4),
)

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# %%
pvals_long = pvals_cleaned.reset_index().melt(
    id_vars=['interaction_name'],  # 或其他标识交互的列
    var_name='cell_pair',
    value_name='pvalue'
)

# %% [markdown]
# ## Fig6: T细胞在TME中的空间定位

# %%
import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import geopandas as gpd
import scanpy as sc
import anndata
import polars as pl
import pickle
import skimage.io
import squidpy as sq

from tifffile import imread, imwrite
from csbdeep.utils import normalize
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point
from scipy import sparse
from matplotlib.colors import ListedColormap
import pyarrow

from csbdeep.io import save_tiff_imagej_compatible
from stardist import export_imagej_rois

# %%
### Here we list BC used for zoom-ins

#                  BC 
#P1CRC  s_008um_00176_00486-1
#P2CRC  s_008um_00460_00073-1 
#P5CRC  s_008um_00365_00224-1 

# %%
adata_merged = adata_merged[~adata_merged.obs['DeconvolutionLabel1'].isna()].copy()

# %%
# 示例使用：
center = "P1CRC-s_008um_00176_00486-1"
crop_coord, bins = get_crop_bins(adata_merged[adata_merged.obs.library_id=='P1CRC'], center, size_microns=500)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")

# %%
is_singlet = adata_merged_merged.obs['DeconvolutionClass'] == 'singlet'
adata_singlet = adata_merged_merged[is_singlet, :].copy()

# %%
sq.pl.spatial_scatter(
    adata_singlet[adata_singlet.obs.library_id == 'P1CRC'],
    color=["DeconvolutionLabel1"],
    library_id="P1CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    shape='square',
    size=0.8, 
    figsize=(8,6),
    img='hires'
)

# %%
# 示例使用：
center = "P2CRC-s_008um_00460_00073-1"
crop_coord, bins = get_crop_bins(adata_merged[adata_merged.obs.library_id=='P2CRC'], center, size_microns=500)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")

# %%
sq.pl.spatial_scatter(
    adata_singlet[adata_singlet.obs.library_id == 'P2CRC'],
    color=["DeconvolutionLabel1"],
    library_id="P2CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    shape='square',
    size=0.8, 
    figsize=(5,4),
    img='hires'
)

# %%
# 示例使用：
center = "P5CRC-s_008um_00365_00224-1"
crop_coord, bins = get_crop_bins(adata_merged[adata_merged.obs.library_id=='P5CRC'], center, size_microns=500)
print(f"Crop coordinates: {crop_coord}")
print(f"Number of spots in crop: {len(bins)}")
sq.pl.spatial_scatter(
    adata_singlet[adata_singlet.obs.library_id == 'P5CRC'],
    color=["DeconvolutionLabel1"],
    library_id="P5CRC",        # 换成你的 sample id
    crop_coord=crop_coord,           # 用刚算好的 crop_coord
    shape='square',
    size=0.8, 
    figsize=(5,4),
    img='hires'
)


# %%
def nuclei_segmentation(img, dim, block_size=4096, min_overlap=128, context=128):
  # Percentile normalization of the image
  # Adjust min_percentile and max_percentile as needed
  min_percentile = 5
  max_percentile = 95
  img = normalize(img, min_percentile, max_percentile)

  # Predict cell nuclei using the normalized image
  # Adjust nms_thresh and prob_thresh as needed
  model = StarDist2D.from_pretrained('2D_versatile_he')
  labels, polys = model.predict_instances_big(img, axes=dim, block_size=block_size,
                                            prob_thresh=0.01,nms_thresh=0.001,
                                            min_overlap=min_overlap, context=context,
                                            normalizer=None, n_tiles=(4,4,1))
  return(labels, polys)


def assign_bin_to_segmented_nuclei(gdf, adata, df_tissue_positions):
  #Set the index of the dataframe to the barcodes
  df_tissue_positions = df_tissue_positions.set_index('barcode')
  # Create an index in the dataframe to check joins
  df_tissue_positions['index']=df_tissue_positions.index

  print("... Creating bin coordinate GeoDataFrame ...")
  # Create a GeoDataFrame for bin coordinates
  geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'],
                                    df_tissue_positions['pxl_row_in_fullres'])]
  gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)

  print("... Overlapping bins with nuclei ...")
  # Perform a spatial join to check which coordinates are in a polygon
  result_spatial_join = gpd.sjoin(gdf_coordinates, gdf,
                                how='left', predicate='within')

  print("... Filtering bins ...")
  # Identify the bins that are in a nuclues and find bins that are in more than one nuclues
  result_spatial_join['is_within_polygon'] = ~result_spatial_join['index_right'].isna()
  barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
  result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlaping_polygons)
  # Remove barcodes in overlapping nuclei
  barcodes_in_one_polygon = result_spatial_join[result_spatial_join['is_within_polygon'] &
                                              result_spatial_join['is_not_in_an_polygon_overlap']]
  # The AnnData object is filtered to only contain the barcodes that are in non-overlapping polygon regions
  filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
  filtered_adata = adata[filtered_obs_mask,:]

  print("... Counting gene expression from overlapped bins for each nucleus ...")
  # Add the results of the point spatial join to the Anndata object
  filtered_adata.obs =  pd.merge(filtered_adata.obs,
                                 barcodes_in_one_polygon[['index','geometry','id',
                                                        'is_within_polygon',
                                                        'is_not_in_an_polygon_overlap']],
                                 left_index=True, right_index=True)

  # Group the data by unique nucleous IDs
  groupby_object = filtered_adata.obs.groupby(['id'], observed=True)

  # Extract the gene expression counts from the AnnData object
  counts = filtered_adata.X

  # Obtain the number of unique nuclei and the number of genes in the expression data
  N_groups = groupby_object.ngroups
  N_genes = counts.shape[1]

  # Initialize a sparse matrix to store the summed gene counts for each nucleus
  summed_counts = sparse.lil_matrix((N_groups, N_genes))

  # Lists to store the IDs of polygons and the current row index
  polygon_id = []
  row = 0

  # Iterate over each unique polygon to calculate the sum of gene counts.
  for polygons, idx_ in groupby_object.indices.items():
      summed_counts[row] = counts[idx_].sum(0)
      row += 1
      polygon_id.append(polygons)

  # Create and AnnData object from the summed count matrix
  summed_counts = summed_counts.tocsr()
  grouped_filtered_adata = anndata.AnnData(X=summed_counts,
                                           obs=pd.DataFrame(polygon_id,
                                                            columns=['id'],
                                                            index=polygon_id),
                                           var=filtered_adata.var)
  return(grouped_filtered_adata)


def plot_nuclei_area(gdf,area_cut_off):
    fig, axs = plt.subplots(1, 2, figsize=(15, 4))
    # Plot the histograms
    axs[0].hist(gdf['area'], bins=50, edgecolor='black')
    axs[0].set_title('Nuclei Area')

    axs[1].hist(gdf[gdf['area'] < area_cut_off]['area'], bins=50, edgecolor='black')
    axs[1].set_title('Nuclei Area Filtered:'+str(area_cut_off))

    plt.tight_layout()
    plt.show()


# %%
img = skimage.io.imread('/data1/home/chenao/reproduction/GSM8594567/spatial/tissue_hires_image.png')

# %%
plt.figure(figsize=(10, 10))  # 设置图形大小
# 使用 plt.imshow() 显示图像
# 如果是灰度图 (2D 数组)，imshow 会自动处理
# 如果是彩色图 (3D 数组)，imshow 也会自动处理
plt.imshow(img)

# 隐藏坐标轴刻度，让图像更清晰 (可选)
plt.axis('on')

# 添加标题 (可选)
plt.title("CytAssist Tissue Image")

# 显示图像
plt.show()

# %% [markdown]
# ## Fig7: Xenium In Situ证实了巨噬细胞亚型和克隆性扩增的T细胞在TME中的存在和定位

# %%
import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns
plt.rcdefaults()
import scanpy as sc
import squidpy as sq
import spatialdata_plot

# %% [markdown]
# 空间数据 io 中提供了一个 Xenium 数据读取器。我们使用它来解析并将 Xenium 数据集转换为 Zarr
# 我们指定数据集的路径以及存储 Zarr 文件的位置

# %%
xenium_path = "/home/chenao/jupyter_home/downstream/process/reproduction/Xenium_data/Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs"

# %%
sdata3 = xenium(xenium_path,n_jobs=10)

# %%
sdata3

# %%
import scanpy as sc
adata3 = sdata3.tables["table"]
adata3

# %%
sc.pp.normalize_total(adata3)
sc.pp.log1p(adata3)
sc.pp.highly_variable_genes(adata3)
adata3.var.sort_values("means")

# %%
sc.pp.pca(adata3)
sc.pp.neighbors(adata3)

# %%
sc.tl.leiden(adata3,resolution = 0.5)

# %%
adata3.obs['leiden'].value_counts()

# %%
sc.tl.rank_genes_groups(adata3, groupby='leiden', method="wilcoxon")

# %%
# 按差异基因进行排序
rg_df = sc.get.rank_genes_groups_df(adata3, group=None) \
          .sort_values(["group","pvals_adj","scores"], ascending=[True, True, False])
topN = 5
top_df = rg_df.groupby("group", as_index=False, sort=False).head(topN)
ordered_genes = top_df.groupby("group")["names"].apply(list).explode().drop_duplicates().tolist()

print(ordered_genes)

# %%
sc.pl.dotplot(adata3, var_names=ordered_genes, groupby="leiden", standard_scale="var", dendrogram=False,cmap="bwr")

# %%
tg = ['Clonotype1_TRA1','Clonotype1_TRA2','Clonotype1_TRA2','STAB1','CXCL10','CXCL9','CXCL11','GZMA','GZMK','PRF1']

# %%
sc.pl.dotplot(adata3, var_names=tg, groupby="leiden", standard_scale="var", dendrogram=False,cmap="bwr")

# %%
cluster_to_cell_type = {
    '0': 'Tumor',
    '1': 'T cell',
    '2': 'Tumor',
    '3': 'CAF',
    '4': 'Macrophage',
    '5': 'Tumor',
    '6': 'Endothelial',
    '7': 'B cell',
    '10': 'Tumor',
    '12': 'Tumor',
    '13': 'T cell',
    # 其他未标记的聚类，统一设置为 'Other/Unannotated'
}
adata3.obs['celltype'] = adata3.obs['leiden'].astype(str).map(cluster_to_cell_type)
adata3.obs['celltype'] = adata3.obs['celltype'].fillna('Other/Unannotated')
adata3.obs['celltype'] = adata3.obs['celltype'].astype('category')

# %%
sdata3

# %%
sdata3.pl.render_images("morphology_focus", channel='DAPI').pl.show(
    title="Morphology Image with 1000-Pixel Grid", 
)

# %%
IMAGE_MAX_X = 30000 
IMAGE_MAX_Y = 32000 
PIXEL_INTERVAL = 1000 # 您希望的间隔

# 1. 获取 Matplotlib Axes 对象
# 确保在 .pl.show() 中使用 return_ax=True 来捕获 Axes
ax = sdata3.pl.render_images("morphology_focus", channel='DAPI').pl.show(
    title="Morphology Image with 1000-Pixel Grid", 
    return_ax=True 
)

# 2. 生成辅助线的位置列表
# 从 0 开始，到图像最大尺寸之外，以确保所有区域都被覆盖
# np.arange 自动生成 [2000, 4000, 6000, ..., 30000, 32000] 这样的坐标
x_line_positions = np.arange(PIXEL_INTERVAL, IMAGE_MAX_X + PIXEL_INTERVAL, PIXEL_INTERVAL)
y_line_positions = np.arange(PIXEL_INTERVAL, IMAGE_MAX_Y + PIXEL_INTERVAL, PIXEL_INTERVAL)

# 3. 循环添加垂直线 (X 轴)
for x_coord in x_line_positions:
    # 只需要确保线在 X 轴的实际数据范围内 (0 到 30000)
    ax.axvline(x=x_coord, 
               color='cyan',  # 换个醒目的颜色
               linestyle='--', 
               linewidth=1.0, 
               alpha=0.6) 

# 4. 循环添加水平线 (Y 轴)
for y_coord in y_line_positions:
    # 只需要确保线在 Y 轴的实际数据范围内 (0 到 32000)
    ax.axhline(y=y_coord, 
               color='cyan', 
               linestyle='--', 
               linewidth=1.0, 
               alpha=0.6)

# 5. 显示最终的图
plt.show()

# %%
from spatialdata import bounding_box_query

def crop0(x):
    return bounding_box_query(
        x,
        min_coordinate=[18_000, 18_500],
        max_coordinate=[19_000, 19_500],
        axes=("x", "y"),
        target_coordinate_system="global",
    )

sdata33 = crop0(sdata3)

# %%
sq.pl.spatial_scatter(
    adata3,
    library_id="spatial",
    shape=None,
    color=[
        "REG1A","SPP1","TGFBI"
    ],
    cmap='inferno',
    img_channel="DAPI",
    img_alpha=0.7,
    wspace=0.4,
)

plt.savefig('/home/chenao/jupyter_home/downstream/process/reproduction/fig/fig7/spatial_expression_xenium.pdf', dpi=300, bbox_inches='tight')

# %%
sc.tl.rank_genes_groups(adata3, groupby="leiden_clusters", method="wilcoxon")

# %%
import matplotlib.pyplot as plt
sdata1.pl.render_images("morphology_focus",channel='DAPI').pl.show(title="Morphology image")

# %%
sdata1.pl.render_images("morphology_focus",channel='DAPI').pl.render_points(
    "transcripts",
    color="feature_name",
    groups='SPP1',
    palette="orange",
).pl.show(title=f"{gene_name} expression over H&E image", coordinate_systems="global", figsize=(10, 5))

# %%
IMAGE_MAX_X = 30000 
IMAGE_MAX_Y = 32000 
PIXEL_INTERVAL = 1000 # 您希望的间隔

# 1. 获取 Matplotlib Axes 对象
# 确保在 .pl.show() 中使用 return_ax=True 来捕获 Axes
ax = sdata1.pl.render_images("morphology_focus", channel='DAPI').pl.show(
    title="Morphology Image with 1000-Pixel Grid", 
    return_ax=True 
)

# 2. 生成辅助线的位置列表
# 从 0 开始，到图像最大尺寸之外，以确保所有区域都被覆盖
# np.arange 自动生成 [2000, 4000, 6000, ..., 30000, 32000] 这样的坐标
x_line_positions = np.arange(PIXEL_INTERVAL, IMAGE_MAX_X + PIXEL_INTERVAL, PIXEL_INTERVAL)
y_line_positions = np.arange(PIXEL_INTERVAL, IMAGE_MAX_Y + PIXEL_INTERVAL, PIXEL_INTERVAL)

# 3. 循环添加垂直线 (X 轴)
for x_coord in x_line_positions:
    # 只需要确保线在 X 轴的实际数据范围内 (0 到 30000)
    ax.axvline(x=x_coord, 
               color='cyan',  # 换个醒目的颜色
               linestyle='--', 
               linewidth=1.0, 
               alpha=0.6) 

# 4. 循环添加水平线 (Y 轴)
for y_coord in y_line_positions:
    # 只需要确保线在 Y 轴的实际数据范围内 (0 到 32000)
    ax.axhline(y=y_coord, 
               color='cyan', 
               linestyle='--', 
               linewidth=1.0, 
               alpha=0.6)

# 5. 显示最终的图
plt.show()

# %%
ax = sdata1.pl.render_images("morphology_focus", channel='DAPI').pl.show(
    title="Morphology image with Auxiliary Lines", 
    return_ax=True # 👈 关键：返回 Axes 对象
)

# %%
from spatialdata import bounding_box_query

def crop0(x):
    return bounding_box_query(
        x,
        min_coordinate=[4_000, 22_000],
        max_coordinate=[6_000, 24_000],
        axes=("x", "y"),
        target_coordinate_system="global",
    )

sdatas = crop0(sdata1)

# %%
fig, axes = plt.subplots(1, 2, figsize=(8, 5))
sdatas.pl.render_images("morphology_focus", channel='DAPI').pl.show(
    ax=axes[0], title="Morphology image", coordinate_systems="global"
)
sdatas.pl.render_labels("cell_labels").pl.show(ax=axes[1], title="Cell labels", coordinate_systems="global")
plt.tight_layout()

# %%
# 让我们可视化这个数据集中用于计算细胞分割的 4 个通道。这特定于使用细胞分割工具包分析并使用 Xenium Analyzer >= 2.0.0 处理的数据集。
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

sdatas.pl.render_images("morphology_focus", channel='DAPI').pl.show(ax=axes[0, 0],title="DAPI", figsize=(3, 3))
sdatas.pl.render_images("morphology_focus", channel='AlphaSMA/Vimentin').pl.show(
    ax=axes[0, 1], title="AlphaSMA/Vimentin (interior - protein)", figsize=(3, 3)
)
sdatas.pl.render_images("morphology_focus", channel='ATP1A1/CD45/E-Cadherin').pl.show(
    ax=axes[1, 0], title="ATP1A1/CD45/E-Cadherin (boundary)", figsize=(3, 3)
)
sdatas.pl.render_images("morphology_focus", channel='18S').pl.show(ax=axes[1, 1],title="18S (interior - RNA)", figsize=(3, 3))
plt.tight_layout() 
plt.show()

# %%
adatas1 = sdatas.tables['table']

# %%
adatas1

# %% [markdown]
# 对于分析，我们使用 anddata.AnnData 对象，其中包含计数矩阵、细胞和基因注释。它存储在 spacedata.tables 层中

# %%
sdata.pl.render_images("morphology_focus",channel='DAPI').pl.show(title="Morphology image")

# %%
from spatialdata import bounding_box_query

def crop0(x):
    return bounding_box_query(
        x,
        min_coordinate=[24_000, 14_500],
        max_coordinate=[26_000, 16_500],
        axes=("x", "y"),
        target_coordinate_system="global",
    )

sdata1 = crop0(sdata)

# %%
sdata1.pl.render_images("morphology_focus",channel='DAPI').pl.show(title="Morphology image")

# %%
adata.obs['region'] = 'cell_boundaries'

# %%
adata.obs['region'] = adata.obs['region'].astype('category')

# %%
instance_id_col_name = 'cell_id' 
region_key_col_name = 'region'

# 覆盖关联，将其指向 'cell_boundaries'
sdata1.set_table_annotates_spatialelement(
    table_name='table',
    region='cell_boundaries',      # 明确指向新的目标元素
    region_key=region_key_col_name,
    instance_key=instance_id_col_name
)

# %%
cluster_to_cell_type = {
    # Tumor 细胞
    '1': 'Tumor',
    '5': 'Tumor',
    
    # Macrophage 细胞
    '25': 'Macrophage',
    '23': 'Macrophage',
    
    # CAF 细胞
    '42': 'CAF',
    
    # 其他未标记的聚类，统一设置为 'Other/Unannotated'
}

# %%
sdata1 = crop0(sdata)

# %%
sdata1

# %%
adata1 = sdata1.tables["table"]
sq.pl.spatial_scatter(
    adata1,
    library_id="spatial",
    shape=None,
    color=[
        "REG1A","SPP1","TGFBI"
    ],
    cmap='inferno',
    wspace=0.4,
)

# %%
xenium_path2 = "/home/chenao/jupyter_home/downstream/process/reproduction/Xenium_data/Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs"
sdata2 = xenium(xenium_path2)
adata2 = sdata2.tables["table"]
adata2

# %%
adata2.layers["counts"] = adata2.X.copy()
sc.pp.normalize_total(adata2, inplace=True)
sc.pp.log1p(adata2)
sq.pl.spatial_scatter(
    adata2,
    library_id="spatial",
    shape=None,
    color=[
        "REG1A","SPP1","TGFBI"
    ],
    cmap='inferno',
    wspace=0.4,
)

# %%
xenium_path3 = "/home/chenao/jupyter_home/downstream/process/reproduction/Xenium_data/Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs"
sdata3 = xenium(xenium_path3)
adata3 = sdata3.tables["table"]
adata3

# %%
adata3.layers["counts"] = adata3.X.copy()
sc.pp.normalize_total(adata3, inplace=True)
sc.pp.log1p(adata3)
sq.pl.spatial_scatter(
    adata3,
    library_id="spatial",
    shape=None,
    color=[
        "REG1A","SPP1","TGFBI"
    ],
    cmap='inferno',
    wspace=0.4,
)

# %%
sq.pl.spatial_scatter(
    adata3,
    library_id="spatial",
    shape=None,
    color=[
        "transcripts"
    ],
    cmap='inferno',
    wspace=0.4,
)

# %% [markdown]
# # 环境信息

# %%
import subprocess
subprocess.run(['conda', 'list'])

# %%
subprocess.run(['pip', 'list'])

# %% [markdown]
# # 单细胞教程全收录
# [B站视频](https://www.bilibili.com/video/BV1S44y1b76Z/):先看一遍视频再去看推送操作，建议至少**看三遍**  
# 以下的资料都可以在这里订阅：
# [scRNA-Seq学习手册Seurat V4修订版](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247510380&idx=1&sn=0eb2b0c519ef1dc25ee8821b736a0fde&chksm=9b3cd03cac4b592ad93149fb4847e6ba3dd03ab69e1ca594ddc64d323fb7496775789cae245d&scene=21#wechat_redirect)    
# [scRNA-Seq学习手册Seurat V5更新版](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247519711&idx=1&sn=d25d4f5daf609bda8367f16db676ecf3&chksm=9b3cf48fac4b7d99e65eba6d4c66c63e40b1d187fa1d9e6449d6cd07d5dec5127cff54886f93&scene=21#wechat_redirect)    
# [scRNA-Seq学习手册Python版](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504877&idx=2&sn=2b56a9cfa19fcc918ee6c2ed88fed017&chksm=9b3caabdac4b23ab6e28425a4a5eb4d40ee3612e51f84f5d1e53ac91bef93f717e9e754398e3&scene=21#wechat_redirect)    
#  
# ## 单细胞测序基础数据分析**保姆级教程**，代码部分整理在往期推送之中：
# [手把手教你做单细胞测序数据分析（一）——绪论](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499776&idx=2&sn=7fc250d94ed88cfd5e2224f071bb1742&chksm=9b3cb950ac4b30467a082527547b91aec15b7d1f40d44e362894aedeb6a99e014518fa5e4743&token=767842863&lang=zh_CN#rd)  
# [手把手教你做单细胞测序数据分析（二）——各类数据结构与读取方法](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499795&idx=2&sn=68e6298eb43633c71939c146d90230ff&chksm=9b3cb943ac4b3055e09739d77e037d679f3d9ec68780d1aa8c592f806fcc8d2322eaee78752d&scene=178&cur_album_id=2923694853375688706#rd)  
# [手把手教你做单细胞测序数据分析（三）——单样本分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500010&idx=2&sn=2f9fd1b0101d579bcd250fdc24e98bfc&chksm=9b3cb9baac4b30ac6f0ce28a6dd10470e73d1e1d6defece5636939a32a6db9b2cce0ff8688eb&scene=178&cur_album_id=2923694853375688706#rd)  
# [手把手教你做单细胞测序数据分析（四）——多样本整合](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500117&idx=2&sn=cef0b0abec6d73a8c2615980ef68aec1&chksm=9b3cb805ac4b3113fdc820b9c9667f9554b2e39f90c45b5806021a415698b0a5141c521c64dc&scene=178&cur_album_id=2923694853375688706#rd)  
# [手把手教你做单细胞测序数据分析（五）——细胞类型注释](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500117&idx=2&sn=cef0b0abec6d73a8c2615980ef68aec1&chksm=9b3cb805ac4b3113fdc820b9c9667f9554b2e39f90c45b5806021a415698b0a5141c521c64dc&scene=178&cur_album_id=2923694853375688706#rd)  
# [手把手教你做单细胞测序数据分析（六）——组间差异分析及可视化](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500117&idx=2&sn=cef0b0abec6d73a8c2615980ef68aec1&chksm=9b3cb805ac4b3113fdc820b9c9667f9554b2e39f90c45b5806021a415698b0a5141c521c64dc&scene=178&cur_album_id=2923694853375688706#rd)  
# [手把手教你做单细胞测序数据分析（七）——基因集富集分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500357&idx=2&sn=40034079b771c478b482e849ed4b216c&chksm=9b3cbb15ac4b3203d517ff71073ab0520fe144bc108b7a26f6b17e88fb63a78c5777946618e3&scene=178&cur_album_id=2923694853375688706#rd)  
# [Seurat中分类变量处理技巧](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486134&idx=1&sn=e581b8ac20bbc968aeef3f9cf56f3a0f&chksm=9b3f73e6ac48faf0f1c81752b0190abc617eb720a71806950e7851c7d19b7e3e092aa97bd005&scene=21#wechat_redirect)    
# [Seurat中分类变量处理技巧](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486134&idx=1&sn=e581b8ac20bbc968aeef3f9cf56f3a0f&chksm=9b3f73e6ac48faf0f1c81752b0190abc617eb720a71806950e7851c7d19b7e3e092aa97bd005&scene=21#wechat_redirect)
#   
# ## 单细胞图片改造计划
# [改造单细胞降维图| 1.DimPlot的探索](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247496162&idx=1&sn=598c2920f7cc0697ee5a1fa6cadaee04&chksm=9b3c88b2ac4b01a415feeebbbec6c4891319439dc8dcdafd33cca7a1222e5adc723ae9d3a73d&scene=21#wechat_redirect)  
# [改造单细胞降维图| 2.ggplot2中的DIY](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247496239&idx=1&sn=1fe6fd4a7fb88ecc0456b4f5710cb5d6&chksm=9b3c8b7fac4b0269e0abeab951ad6ae0ddd1464dd9d1867280eeb9755fda5ba8f5237d88a7a6&scene=21#wechat_redirect)  
# [改造单细胞降维图| 3.3D降维图与动图绘制](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498080&idx=2&sn=a9bcf65e3ae34ba201e8d49ee33c0728&chksm=9b3c8030ac4b0926c3f059a9a5d23f1bcdbcaf8e5a562eba2b7e3bec3396c44573269bd08e7b&scene=21&cur_album_id=2707677276481814529#wechat_redirect)  
# [改造单细胞降维图| 4.一些造好的轮子](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498191&idx=1&sn=952312a6b26232d460bf787cc276543f&chksm=9b3c809fac4b098972c8c1dd83b9f954a38fe5b5f6bc7e295f3f3a7ba6c4c3d32951679ac24d&scene=21#wechat_redirect)  
# [沉浸式统计细胞比例](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503971&idx=2&sn=7dafb1e743ebbb95e488d8d1b7024448&chksm=9b3ca933ac4b2025f922a325ec6893e93613f130ff1d4fd609e2d7274f7f1fdb5134092e28ad&scene=21#wechat_redirect)  
# [改造单细胞FeaturePlot](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504616&idx=1&sn=f7b4a0e6b25d224380cd803e3463a3e4&chksm=9b3cabb8ac4b22aeb65e0d6ec32a1f8674dcf907da078ccb3b61e06c4aa035864ebec092f571&scene=21#wechat_redirect)  
# [VlnPlot用法及ggplot复现](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504768&idx=1&sn=4fc5eb47037d3a80178e32975488d70b&chksm=9b3caad0ac4b23c61147cfb0f46ecd846568f37ce89c963eb5f069d3e89a8ba285897726bb3d&scene=21#wechat_redirect)  
# [VlnPlot的ggplot改造](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504843&idx=1&sn=b795a3e36ae092d8f87abb5a87187828&chksm=9b3caa9bac4b238dd33e905c33924b58aa4ff14e83abbf08046fba7fe765ecb0204163cb6acf&scene=21#wechat_redirect)  
# [改造单细胞DotPlot](https://mp.weixin.qq.com/s/GPp-QYwT6jOjFfeZjs9hZg)  
# [改造单细胞DoHeatmap](https://mp.weixin.qq.com/s/V_JhpBOy1wG6U-uQsj8CPQ)   
#
# ## 3、SCENIC转录因子分析：
# [1.绪论](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500098&idx=1&sn=8c31241bd259e303a6ac3c9740f298d9&chksm=9b3cb812ac4b31046fd298c8ee5f0bc0de445df62b234bd934c34ca655e61de3da54d81e2053&scene=21#wechat_redirect)  
# [2.学习手册](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500439&idx=1&sn=d6b4874a8352f5b7a4165629c03058ce&chksm=9b3cbbc7ac4b32d158fd3c5c2f155c54c892d40eaae3718409936fc319f69705da32badf79d6&scene=21#wechat_redirect)  
# [3.软件安装与数据准备](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500631&idx=1&sn=5d73db4dcc74b66de07dee9db910bb83&chksm=9b3cba07ac4b3311d1820bb305778aea707294f19fd6ff2631cd003da80bdcff8f8bc8ed34de&scene=21#wechat_redirect)  
# [4.精简版流程](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501260&idx=1&sn=f1e2a165599e753e2959aad7c8a298cb&chksm=9b3cbc9cac4b358a9499654a2f95fa9010309ee4f408ef8b1f31ceb1a0f263ce22e8b511c430&scene=21#wechat_redirect)  
# [5.step1+step2构建共表达网络与regulon](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501692&idx=1&sn=3f7e263a93942ad03faa1ba5f8694ec4&chksm=9b3cbe2cac4b373af9db74a50ffb08a9b8a3b958fc6385c937385f0f7323bc935cdb0fe5275c&scene=21#wechat_redirect)  
# [6.Step3 利用AUCell对Regulon评分](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501989&idx=1&sn=2700bf779d7d261a4d751d1635fdce2c&chksm=9b3cb1f5ac4b38e3014476242fe956cf510cc45b9f4cef86070af8094fd5f884f9d5c3a3363f&scene=21#wechat_redirect)  
# [7.Step4 二元矩阵的计算与可视化](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503995&idx=1&sn=26c7e9ec56af434e62b02853cfe667d6&chksm=9b3ca92bac4b203d7986645ab5a3fb574f44185c7ced8bd47957fa237d436a95bfa2579d052a&scene=21#wechat_redirect)  
# [8.Step5 regulon聚类、分群、降维](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504024&idx=1&sn=143c0ec4d94b8fca0a4a47d3dc24415b&chksm=9b3ca9c8ac4b20de70a87acb2d3f810b0e1497071b1c6f14f1882a6d844fff016e9ba5581eb3&scene=21#wechat_redirect)  
# [9.下游探索](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247505276&idx=1&sn=949c753feb85d5877121433caaa33013&chksm=9b3cac2cac4b253affb587a937d3f5e0576851a4f4e1570cb7510a594d155873961bbd91758d&scene=21#wechat_redirect)  
# [SCENIC转录因子调控网络图](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247520750&idx=1&sn=496bead1531a7eddd5c286030d350745&chksm=9b3ce8beac4b61a8589a1f8a66bd0571e09c0646655b6d805348bf6f4a0c7b52851615c22868&scene=21#wechat_redirect)  
#   
# ## 上游**fastq**文件处理：
# [单细胞分析的最上游——处理Fastq文件：cellranger](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484923&idx=1&sn=b5876af14fbee68d1e9db4b0f70cd1c8&chksm=9b3f7cabac48f5bdb3aac7d20201d89121a83720b41a10017564fedbb8c4a13df0da92daffd7&scene=21#wechat_redirect)  
# [单细胞分析的最上游——处理Fastq文件：dropseqRunner](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484937&idx=1&sn=4fe1fa518202bdaf462bb93a2614ed32&chksm=9b3f7f59ac48f64fafaf34bac2c80f06d3a86b6a84ccf12566c7eac25bb71c7edd3d74004cf7&scene=21#wechat_redirect)    
# [Cellranger报错：Unable to distinguish between](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504739&idx=2&sn=e4e6ad0559b9b6a703d6ab37c2a8fdb0&chksm=9b3caa33ac4b23257424184bda44dd600ad98bc0dfdfa90ab0e3ebaf208e598b0b710fbafa4b&scene=21#wechat_redirect)  
# [BD Rhapsody平台单细胞转录组定量流程](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506542&idx=1&sn=99cf7af7fc0d56b0d6fe7ea3596df2b6&chksm=9b3ca33eac4b2a28e6682a7cf0b2091a387a06ecf89ee10737cf2194286ec0bba472ef1c2b90&scene=21#wechat_redirect)  
# [BD Rhapsody定量输出文件读取](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247519706&idx=1&sn=9d084bd8edf6ab4d9184672e5a4a2722&chksm=9b3cf48aac4b7d9c16ef85de7a5f33e65552bbba5d11d2170237063cf754fe9ff376647d3bd7&scene=21#wechat_redirect)  
#   
# ## **细胞通讯**
# [跳转B站，连播更方便](https://www.bilibili.com/video/BV1Ab4y1W7qx?p=1)  
# **往期推送**  
# [《细胞通讯》CellChat学习手册](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501735&idx=2&sn=35fc792961507f4ae1cbae3a2882fba6&chksm=9b3cbef7ac4b37e1417d4108c8bd378c68c3907431ddfc52ad37c916e85649bb450d761a8e19&scene=21#wechat_redirect)  
# [《细胞通讯》1.概论](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501282&idx=2&sn=cd37ea80e0b92f5bfbc98ea286e363da&chksm=9b3cbcb2ac4b35a4258d21fe003374795619fc2356376ebbe4ff640233f4b4c6e3e3ead41370&scene=21#wechat_redirect)  
# [《细胞通讯》2.1CellChat基础分析教程](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501296&idx=2&sn=fb3e252d95ee8b50cd446e71335b1aa0&chksm=9b3cbca0ac4b35b621b4f51d2c3a130f00461faeeff067396222c2fc2eeb53f3544664d7d9bc&scene=21#wechat_redirect)  
# [《细胞通讯》2.2CellChat多组别分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501523&idx=2&sn=6fa803f256c2e0b945f93ae7d8245f50&chksm=9b3cbf83ac4b3695b587e9817840f419104c8680f5ef08ee39aa149b79852f48b179f68e3397&scene=21#wechat_redirect)    
# [CellChat空转细胞通讯合辑](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247512011&idx=1&sn=fbd41832c5995d79fead016931d3bfbd&chksm=9b3cd69bac4b5f8d16e7d8b8eeb5490e4938daf68902b06a8d60e25f150ec4ec459796f3aaac&scene=21#wechat_redirect)  
# [SeekSpace细胞通讯分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511890&idx=1&sn=3457b875cd00b4433f50a39a45e1cb8a&chksm=9b3cd602ac4b5f144ad982dfbdf583b56d1bfb768e4e5f2cc9e5424af6e014e752835fadc8a6&scene=21#wechat_redirect)  
#                               
# ## **拟时序分析**
# [B站连续播放起来比较方便](https://www.bilibili.com/video/BV1br4y1x7Hf?p=1)  
# **往期推送**  
# [1.概论](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486553&idx=1&sn=b689ec962a17ca1a1ada29ae8ff8c019&chksm=9b3f7509ac48fc1f0f81e288de8f15c86548215a1cd9d5c740dc63f2421371d13890398ab1df&scene=21#wechat_redirect)  
# [2.monocle概论](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247487001&idx=1&sn=329129e75f703e3f332e76af9c66871e&chksm=9b3f7749ac48fe5f081248c2a48c56685e3915015036c73f8f67ec3e6acc246ddf815daedb22&scene=21#wechat_redirect)  
# [3.monocle2实操：精简版](https://mp.weixin.qq.com/s/k-1Jdc0J9EXb2jfwo7935w)    
# [4.monocle2实操：完整版](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497730&idx=2&sn=930dc0e276219d70c937e943119388ff&chksm=9b3c8152ac4b08446d457cef3ede4a4a17c0ebacaec0c36efe7d7c8706fb4a3d62222d702ee5&scene=178&cur_album_id=2820794880947159043#rd)  
# [5.初识monocle3](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247489044&idx=1&sn=061dbfc34d200b73a8227d0bb28a1d27&chksm=9b3f6f44ac48e65251f7e5ecce84637c05d50ecb92769e8386f221f8e7b9bfeb719232def863&scene=21#wechat_redirect)  
# [6.monocle3的降维、分群、聚类](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247489494&idx=1&sn=44f0a44857c12fcf9a1ab84e4f2cbebf&chksm=9b3f6e86ac48e79026321de2fdca81a139d104d655884ffe87a68e0a1551ddffa5af7919d19e&scene=21#wechat_redirect)  
# [7.monocle3的拟时序分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490418&idx=1&sn=18ae0b6139fe1a2436dcf6cc47850297&chksm=9b3f6222ac48eb3450de0b02a34eb5fff2db89ab8f155460d554cbb23be0cd27fd349079e7b2&scene=21#wechat_redirect)  
# [解决monocle2的orderCells报错的两种方法](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491083&idx=1&sn=c5697ab024193cb6595d27e4f4d3bfaf&chksm=9b3f675bac48ee4d7244d9819e14fe210f7e42e40fc7a6a8df7b95bfd92f06334d60d226cb93&token=333450324&lang=zh_CN#rd)  
# [一文搞定拟时序分析的下游可视化探索](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497634&idx=1&sn=64872660e8ce80ff08c6c0745c2b241b&chksm=9b3c8ef2ac4b07e4a2ef03870c91fdc1c3ded245dcecd5a2049841cb7e9749bd1cd0852133f1&scene=178&cur_album_id=2820794880947159043#rd)  
#   
# ## 单细胞拷贝数变异
# [单细胞拷贝数变异学习手册](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522733&idx=1&sn=4e3a0707fd8f83b5ad9377bc91ee8e53&scene=21#wechat_redirect)
#
# ## scATAC
# [scATAC分析学习手册](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247521636&idx=1&sn=208907010f3f7412c1c91343d7dfb0d3&scene=21#wechat_redirect)  
# [一文了解单细胞染色质及常见分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522998&idx=1&sn=843d378d3b679951ad3c8e5da7c32d87&scene=21#wechat_redirect)  
# [一文了解染色质可及性和表观基因组调控](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247510114&idx=1&sn=f7e8b1df1804a6bba0d79b1ffec2467c&scene=21#wechat_redirect)  
# [snATAC联合scRNA揭示小鼠肾脏细胞分化程序和疾病靶点](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522324&idx=1&sn=62993c2354af1e793fcfd0e06586af17&scene=21#wechat_redirect)  
# [百万级单细胞多组学数据集成](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499336&idx=1&sn=864c535b38ce630487fccc602d893a80&scene=21#wechat_redirect)  
#
# ## **scFAST系列**
# [01.绪论](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503875&idx=1&sn=2b5ca0c57e17e51abb52adfe199c5736&chksm=9b3ca953ac4b20450f643b94bcfc95aa2b04ded197ce93e455a281901d7ccba099b48cfd03d9&scene=21#wechat_redirect)  
# [02.单细胞分析准备工作](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504023&idx=1&sn=4091a22e853242b4aae7a7d3b4008f5c&chksm=9b3ca9c7ac4b20d1bb0389fad5fe85a3c4cf49179747fe97dffaab1087b7d2793668d9f63aa4&scene=21#wechat_redirect)  
# [03.fast模块](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504044&idx=1&sn=d37991280217dcc29280efa344f856bb&chksm=9b3ca9fcac4b20eab124ccd4ba6ad3f5eee640df314e7d8fec0225035f06178e3a51dbe2033f&scene=21#wechat_redirect)  
# [04.Seurat分析及可视化](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504435&idx=1&sn=c6a6f7ca00684601f3cda49e53be5339&chksm=9b3cab63ac4b2275d993bbaa64c71652685065fa377be57319fc6a5c47eb958b3644c7459cd5&scene=21#wechat_redirect)  
# [05.mut模块](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504878&idx=1&sn=2c5d9fc694eaf3e264c38b0691cd3a43&chksm=9b3caabeac4b23a87d201009137572c02ba1ea0dd747bcd8b8c02878b45917b6e235551c41cb&scene=21#wechat_redirect)  
# [一文学会scFAST-seq数据分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247505213&idx=1&sn=5f4ed0587dd96149d127a98d36f721f0&chksm=9b3cac6dac4b257bbc00c4647382d27883139bf06ab3071947bf1cccde8d1d501006dabab5f1&scene=21#wechat_redirect)  
#   
# ## 其他单细胞相关**技术贴**:
# [应审稿人要求| pseudo bulk差异分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247521437&idx=1&sn=d47b48170f7189b8fb74b2316d0c4898&chksm=9b3cedcdac4b64db43487176203e52e53fe08e01ce819a1615fa1e2f9f8d513bbd4efdee5558&scene=21#wechat_redirect)  
# [genesorteR——比FindAllMarkers快一百倍！](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247521266&idx=1&sn=2a74f282222b694f51362445e65c85a6&chksm=9b3ceaa2ac4b63b4a721d1134317af5ee1c9fbe844383a1f404249748e506b674ea26e66d0d2&scene=21#wechat_redirect)  
# [单细胞必修课——Seurat流程](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247512220&idx=1&sn=c71bc3f1ef007eb29be8b3b4ec10f0af&chksm=9b3cc9ccac4b40da749d2388dc2d885239a1e4f6ea0989f36b031c7253ce103b061e2605ad8c&scene=21#wechat_redirect)  
# [scRNA-Seq双细胞过滤手册](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511622&idx=1&sn=60b35668e569267bab14ec8fbba8bd74&chksm=9b3cd716ac4b5e00122eac47f55085198040a4a0e97c56895bbc39f3a9f1a1eb4eba27048a8f&scene=21#wechat_redirect)   
# [scRNA-Seq学习手册Seurat V4修订版](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247510380&idx=1&sn=0eb2b0c519ef1dc25ee8821b736a0fde&chksm=9b3cd03cac4b592ad93149fb4847e6ba3dd03ab69e1ca594ddc64d323fb7496775789cae245d&scene=21#wechat_redirect)  
# [快速上手Seurat V5](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247509880&idx=1&sn=3b7cd169ccfe688d10cb1ca70dc91e91&chksm=9b3cde28ac4b573ec6d5f4bceafa4294c0e9a001b7d0c0b64f274a9978b1438eeccc747ca0c0&scene=21#wechat_redirect)  
# [单细胞Figure 1常见图表(视频教程)](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247509260&idx=1&sn=9ff90a428c3f22d43a31ec864e2d8a94&chksm=9b3cdc5cac4b554a2f7cb2a6368f7df9a9a3f6e406eb070663f53d81841a76e188c0b50727de&scene=21#wechat_redirect)  
# [SCTransform与ScaleData](https://mp.weixin.qq.com/s/6jwfDa22EUgo2wYpmxCCFg)  
# [细胞的数量由誰决定？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484056&idx=1&sn=4967a054f5241ff3a1bd4a6fe9a502b1&chksm=9b3f7bc8ac48f2deca3141461be035deb8658be184d78597af28f8ba2092d12f020768275248&scene=21#wechat_redirect)  
# [单细胞中应该如何做GSVA？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484701&idx=1&sn=ad90680f252d092a2e21f349c8f32cab&chksm=9b3f7c4dac48f55bf1d11337406ab0e1295e8e1731dbac0715eedc99d11beff7f380e42dd49b&scene=21#wechat_redirect)  
# [答读者问（三）：单细胞测序前景](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484280&idx=2&sn=844af02070012a1c2681880596ba51a2&chksm=9b3f7a28ac48f33e95eb906121be8953410ca4c572c8035cb801609a56caf9dcff57fa989602&scene=21#wechat_redirect)  
# [答读者问（四）：如何分析细胞亚群](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247484316&idx=2&sn=e8d02fddc80d1137bb64c55541f1651d&chksm=9b3f7accac48f3da7df7fcad9c996291bc66967c04d16dde92ed391c56e610a7298ac8fb7dbb&scene=21#wechat_redirect)  
# [答读者问（八）：为什么Read10X也会报错?](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485012&idx=1&sn=6007c8d7cb2a236dd7240a678dac0391&chksm=9b3f7f04ac48f612f811c1df9e9450dc8112d343e4bf245be5149c2e2a073d3e339cd0e3ff61&scene=21#wechat_redirect)  
# [答读者问（十）整合后的表达矩阵，如何拆分出分组信息？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485134&idx=1&sn=13f5f8c47d8af2beeca4e103a5f5e26d&chksm=9b3f7f9eac48f688b454a1062ee4e41f7072ac9d58e5466fa7ed4bb05018adbfe167bb3a445a&scene=21#wechat_redirect)  
# [答读者问（十一）如何一次性读取一个目录下的cellranger输出文件？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485145&idx=1&sn=402cb28c16aba0dc647037c9b7aba8f7&chksm=9b3f7f89ac48f69fe27cd427df3eff2e88c240bf20dfe92d8b522eecd7302d132d880bf29047&scene=21#wechat_redirect)  
# [给你安排一个懂生信的工具人（十）：不学编程 零代码完成单细胞测序数据分析：Loupe Browser](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485730&idx=1&sn=42b9b35026bec5359ac34750dea3bf44&chksm=9b3f7072ac48f964970f673a86e27ecd453a80b953168dc110e0ea6bfa4110eeed3350c9601a&scene=21#wechat_redirect)  
# [什么？不做单细胞也能分析细胞类群和免疫浸润？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485798&idx=1&sn=c849bc7a0093af92fc99e13b55842bb2&chksm=9b3f7036ac48f92048d726b38f2dd978d0b0cce779c1a4090aa4da5cb2df3ad30a5c5dbed27e&scene=21#wechat_redirect)  
# [答读者问 （十三）查看Seurat对象时的ERROR：type='text'](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485859&idx=1&sn=a014c4cb6034551e57c70820c7781e49&chksm=9b3f70f3ac48f9e52620e3ff197eb98e58836e8bb6baad634387a521fdda1e4a194356517c73&scene=21#wechat_redirect)  
# [各类单细胞对象（数据格式）转换大全（一）](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485925&idx=1&sn=14e08440200f0f5cb4dcc2dc615efeb8&chksm=9b3f70b5ac48f9a37a43cc2a8b2a5e485a3f84b4ef87949dfe99056bc6301b110e80800f6bb1&scene=21#wechat_redirect)    
# [单细胞对象(数据格式)转换大全|2. h5ad转Seuratobj](https://mp.weixin.qq.com/s/iYz4ZsD-oJ2ZJsZpbGKRsg)   
# [批量整理好GEO中下载的单细胞数据](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486045&idx=1&sn=0b6874bcffc8abafb505cabce4c4fc0e&chksm=9b3f730dac48fa1b140ef01c07a7749e69df4fa7d2056deeed0673c668b6246724cd85741125&scene=21#wechat_redirect)  
# [答读者问 （十四）Seurat中分类变量处理技巧](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486134&idx=1&sn=e581b8ac20bbc968aeef3f9cf56f3a0f&chksm=9b3f73e6ac48faf0f1c81752b0190abc617eb720a71806950e7851c7d19b7e3e092aa97bd005&scene=21#wechat_redirect)  
# [答读者问 （十五）稀疏矩阵转matrix, as.matrix函数是下下策](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247486259&idx=1&sn=97de5c9e15e48cb62c8a680376a4ee0d&chksm=9b3f7263ac48fb75a91b40416826fe35b340284cf38255925ac447443a62cdc1df082d20fb61&scene=21#wechat_redirect)  
# [答读者问 （十六）做单细胞测序到底需要多少内存](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247487923&idx=1&sn=15043c6d729b0d5493ab875b0e54a259&chksm=9b3f68e3ac48e1f5dae514a01c857f9a56a90e2977dd7862db28d06fe55888cd9bb4afbb14e0&scene=21#wechat_redirect)  
# [答读者问 （十七）调用的线程越多就算的越快嘛?](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247488995&idx=1&sn=7140e0d037671542fb37dff3a3c94b28&chksm=9b3f6cb3ac48e5a5f54a02be9e6af9f5a4cf8da1de1770f9035130c84dedf059a024964fc15a&scene=21#wechat_redirect)  
# [答读者问（十八）、一个我至少被问过30遍的monocle报错](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490010&idx=1&sn=24c55fa4bd5d7d829625cb67dbaf321a&chksm=9b3f608aac48e99c9f27e3423c8398034a1819d0ccebbf21170b955461d2c335132d6f3cce25&scene=21#wechat_redirect)  
# [一文搞定单细胞基因集评分](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497831&idx=2&sn=82fa51c76a381fba0b444e9b6b23eb5f&chksm=9b3c8137ac4b0821dfb62e9b365ea3c32abbbd80df42df63181c37ad9ed991c0183c2c7981d8&payreadticket=HCkRKr7R3-O3p5-0VkAwOOfgIXR8STjX9HoZYYWKS3jgJtptsG2Mu4KWI6Pz6cAuCOxiFPM#rd)  
# [没有barcode文件的单细胞数据要怎么读取](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490952&idx=1&sn=364a3941e59e45234795c906a42fc9ab&chksm=9b3f64d8ac48edcec70384cf5c6a91224d0f7841ce3a7cfcade6513500c13673a6695d6f0299&token=1760636281&lang=zh_CN#rd)  
# [单细胞基因集评分之AUCell](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491761&idx=1&sn=41f5c1ceaacf6a1c4840567f8a349e06&chksm=9b3c99e1ac4b10f7ca158561e96fec76a49ff32b8412cb146ed3fb7aa9260c75172b4eb16f9a&scene=21&token=7173752&lang=zh_CN#wechat_redirect)    
# [如何加快Seurat的计算速度](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494116&idx=1&sn=f3dc577cbb464779062aada313191a05&chksm=9b3c90b4ac4b19a27a17f958e3437337d5b10b6ba7153e4dd43dedcb4df812d8c6b4d86fe970&token=416896563&lang=zh_CN#rd)  
# [粉丝来稿|1. Seurat4相较于Seurat3的几点改动](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494516&idx=1&sn=cb2d4e6b182428ce7b98e7fff2e8729a&chksm=9b3c9224ac4b1b32c4fe99758960305c4845395cfe03851d70a0407043ee36250e55d20106fc&scene=21&token=1634975705&lang=zh_CN#wechat_redirect)  
# [如何加快Seurat的计算速度](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494116&idx=1&sn=f3dc577cbb464779062aada313191a05&chksm=9b3c90b4ac4b19a27a17f958e3437337d5b10b6ba7153e4dd43dedcb4df812d8c6b4d86fe970&token=468535359&lang=zh_CN#rd)      
# [答读者问（二十）四个单细胞样本只给了一套文件怎么读](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494879&idx=1&sn=c2287ab3684792747dfac53d6704462b&chksm=9b3c958fac4b1c994a4b075c4f8475816e121966fa8c52cb0a2ce8c14987988ab9bdb49bdf69&token=468535359&lang=zh_CN#rd)
# [人类单细胞测序数据中有哪些以"**-"开头的基因](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247495019&idx=1&sn=f52ac22ca1ea47ecddaa3503bdc2c0b6&chksm=9b3c943bac4b1d2d521f7155ce8f499389dade23828017572a84cdce71e5d08d37d74e08512b&scene=21&token=1624480671&lang=zh_CN#wechat_redirect)  
# [为什么总把分辨率调的很高](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247495765&idx=1&sn=2ec5e3e26491e28c92e25c40630bffae&chksm=9b3c8905ac4b0013a23f9cda2cccd67cdf3c4a941bac161e3a15a12ab98cd4e5eca8aeedbe96&token=133749967&lang=zh_CN#rd)            
# [答读者问（六）：Seurat中如何让细胞听你指挥](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497669&idx=3&sn=ac3408c7dd28425426d954f81d283252&chksm=9b3c8e95ac4b0783c12182a960e8158c8a0992a20335f8c70f1dc63855c5984dd22e4dee706d&scene=178&cur_album_id=2820794880947159043#rd)  
# [答读者问| 22.object 'CsparseMatrix_validate' not found](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247496273&idx=1&sn=fd5a92dd7fdc349b1624784ca74b1220&chksm=9b3c8b01ac4b0217cf06d89bfea49c3fa0c9b580957338c0cc8f55ede3e8f65d6d7162ae27e3&scene=21#wechat_redirect)  
# [单细胞数据在R中的读取及存储速度太慢怎么办？这个神包来帮你！](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500270&idx=1&sn=11bdf3c7f9042fc03f9ed237a9b9677b&chksm=9b3cb8beac4b31a8bf1d9fb862553ed32ac9a56f4a6e0f2d38c9478132f87327d8c6b6940573&scene=21#wechat_redirect)  
# [单细胞分群分辨率选择困难症](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501241&idx=1&sn=9e15d23f721ad7f4ae1335ce5362031b&chksm=9b3cbce9ac4b35ff43434c0384afbbf5689f5b62a677e55e8a635d4812da57593d85a73cf48b&scene=21#wechat_redirect)  
# [一文搞定空间转录组与单细胞测序的整合分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507210&idx=1&sn=f32037ffefd60b0d1bcca7bf9c3e0823&chksm=9b3ca45aac4b2d4cf0c68299db3a2891130d89d9e0fe2a9ecb34b394736550bcf7d7c51090cc&scene=21#wechat_redirect)  
# [一文掌握十个单细胞数据库](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507979&idx=1&sn=e94225b592733f2e22b4b09245027365&chksm=9b3cd95bac4b504df8dc25d9c5cfddf56d3b47e7cc6f33a7ef369bc916c60abb6b49277b79f7&scene=21#wechat_redirect)  
# [如何将已经构建好的Seurat空转对象转存为原始文件？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247508133&idx=1&sn=297b8fb9cc8f33a1d19167dba1104e97&chksm=9b3cd9f5ac4b50e36df53e3e18bb2990d9635d88a725f1e54fc5eb5967ac9f0934b79cbb8b34&scene=21#wechat_redirect)  
# [一文学会零代码单细胞分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511007&idx=1&sn=e7be3d9775890e5717889714629cfaf5&chksm=9b3cd28fac4b5b9936428d5f74d8b581b5cddc18d0fc42e2a5d75604cce080588d50f1304a17&scene=21#wechat_redirect)  
# [SeekSpace| 会单细胞就会空间转录组](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511779&idx=1&sn=23868525d86a9ed7867a3837ed18612c&chksm=9b3cd7b3ac4b5ea502950d0278d68da978f6011b2597c3ec90fb77e79320883fb4a5be2d5399&scene=21#wechat_redirect)  
# [10xGenomics官方可视化工具——Loupe Browser](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247519698&idx=1&sn=0135ef1763fc146252835f2bd7f728ec&chksm=9b3cf482ac4b7d949d0936f11b487d5da4fd8f235bb86f63e960e23785885b9db0d14897b63d&scene=21#wechat_redirect)  
# [快速上手Cellranger](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247521051&idx=1&sn=e2ab820001013bae6e35310e27897086&chksm=9b3cea4bac4b635ddbb74c29240896902b497cdc47f1f37d99fa0fb8a402375e513690f24e3f&scene=21#wechat_redirect)    
# [单细胞/空转生信工程师自测卷(入门级)](https://mp.weixin.qq.com/s/2RvmlET91UkS-13pu3ZrVw)
# [30w单细胞数据会吃掉多少内存？](https://mp.weixin.qq.com/s/0WoO08Yb80g7sa_R8s1ZwQ)
#
#               
# ## 单细胞文献阅读
#  [《Advanced Science》新型Arf1抑制剂促进癌症干细胞衰老并增强抗肿瘤免疫](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522233&idx=1&sn=287a380de741ff55f0ce78dc9c727985&chksm=9b3ceee9ac4b67ffd925abafeba2cb215a5984e6dc85f0d8607340a07ec9208f1d19eae0e04b&scene=21#wechat_redirect)  
# [IF14.3| scRNA-seq+脂质组多组学分析揭示宫内生长受限导致肝损伤的性别差异](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247520900&idx=1&sn=157134292ffc2e5df730c5ba8a357b61&chksm=9b3cebd4ac4b62c286543c30381e9b4fdb8b44b28d4fa902fb963bbd45b28e1d3f1bc98eb148&scene=21#wechat_redirect)  
# [Q1临床设计+生信分析，单细胞揭示骨髓增生异常综合征中的造血功能异常&微环境中的细胞状态](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247520295&idx=1&sn=8f98781c4ecc8205bcfd2cfa86caed9c&chksm=9b3ce977ac4b60617583775dcd4f0b2e3133a4fec23f7c9b95ae42b26bae9b8f6ddcd1f5ccb0&scene=21#wechat_redirect)  
# [13分+文章利用scRNA-Seq揭示地铁细颗粒物引起肺部炎症的分子机制](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247511886&idx=1&sn=31c712d00377b0f0ed33703b5b4fbfd8&chksm=9b3cd61eac4b5f084fde14a2eb50356f2f40a66c608d38909f12a3d81f48cff029a09ce2ee8c&scene=21#wechat_redirect)    
# [除了铁死亡，还有铜死亡？！](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247518144&idx=1&sn=d036f1bb04b09eb73f53e2104d69404c&chksm=9b3cfe90ac4b77862f2f463b42132c7bfdb687a1b61d350dcf656117b1331271e179710fd4fa&scene=21#wechat_redirect)  
# [scRNA-seq揭示脓毒症预后水平预测的关键靶点！](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522281&idx=1&sn=f95eb12fdc340b8da1123ec6253e6fac&chksm=9b3ceeb9ac4b67afa587e5fa3a2a8960e9bb28a1ad6dc48e49120a76154417a4c8be95af90b0&scene=21#wechat_redirect)  
# [《Narure Biotechnology》教你怎么"捏造"单细胞数据](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247521498&idx=1&sn=925561221f12e2ef2846b897914cc692&chksm=9b3ced8aac4b649c9073ec8dcaae49d689acfaf3d2d8b0012b727ef75b8231088b4e176fa9ca&scene=21#wechat_redirect)    
# [拟时序分析神包—monocle的三篇《Nature》](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247519064&idx=1&sn=173948cb696c7d33b90cd21b28986777&chksm=9b3cf208ac4b7b1e328fdc445f9aecd6cfb836fdb33767876ba1aa35c55098a2a071dbe8e2f0&scene=21#wechat_redirect)  
# [测序技术的发展与应用](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247495792&idx=2&sn=1aecc80b68f65f8a44a899df524038c0&chksm=9b3c8920ac4b0036b23268331b457a335b22b4fa3aed1db287ab79e7a4613df8d195bd11421a&scene=21#wechat_redirect)  
# [单细胞测序解析糖尿病肾病中肾小球的动态变化](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497682&idx=2&sn=7ef49191c5e334bb291ccdaa8ea444ae&chksm=9b3c8e82ac4b079425433938c831d528eddb9237efaca49acbc210e9cd715d057541589dc290&scene=178&cur_album_id=2802463431357661184#rd)  
# [单细胞测序技术解析健康人与T2D患者的胰岛差异](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497713&idx=1&sn=b7f8df758bc4d95b88969378c6dfbeee&chksm=9b3c8ea1ac4b07b7c6a0b5023ebf6c8ecd7d34741ab076baa0b564535d1d753bbf9dc5be5290&scene=178&cur_album_id=2802463431357661184#rd)  
# [小鼠全肾单细胞测序开篇之作](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497713&idx=2&sn=554d6a4e1b28b2e5303695fae5942a34&chksm=9b3c8ea1ac4b07b72323ce25cd49474a9228ec1aaa56243f57dca3f2f668dc42ab9cd58dca21&scene=178&cur_album_id=2802463431357661184#rd)  
# [Nature也能白嫖？](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497867&idx=2&sn=caeb8ed06ffd5239b58cebcd90843cf8&chksm=9b3c81dbac4b08cd9eec50978951a1b890ca08edec644ff75205a09ad1f067dff25d6ecb938e&scene=178&cur_album_id=2802463431357661184#rd)  
# [一篇不花钱就能白嫖的文章](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247497840&idx=2&sn=1dacbfb190e11d2ff104e32f238cf8c4&chksm=9b3c8120ac4b08367bfe3a0dd8d55f21f5d13528826de3144a9a85ce9b1c9f4c710de703edae&scene=178&cur_album_id=2802463431357661184#rd)  
# [高氧下小鼠肺发育损伤的ScRNA图谱](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498513&idx=2&sn=79845edf56a3023c59c73403d7b6494b&chksm=9b3c8241ac4b0b57793402686edc8c35aaa3b0959d3096b5122207fbb70ab809d1df5faee0d7&scene=178&cur_album_id=2802463431357661184#rd)  
# [文献阅读（十二）、IgAN & STRT-Seq](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485376&idx=1&sn=02e228d0993d8a653c8bb1d531021cdb&chksm=9b3f7e90ac48f78686d6d5953aeb192db4f28bf577dbebd76d9fabbe797f135b7fb08fb42867&scene=21#wechat_redirect)  
# [老树开新花——EGFR、肿瘤、免疫+scRNA-Seq](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498596&idx=2&sn=8f63f5856cdf90e15f1d19b4cfe6a39e&chksm=9b3c8234ac4b0b221a5d2c417558f0e63f742459c0a3f8c1cbb4cdf74ce25ddcaf0e48f4a7ce&scene=178&cur_album_id=2802463431357661184#rd)  
# [癌前基质细胞驱动BRCA1肿瘤发生](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498597&idx=2&sn=9f5fee38d0821052c9c6b3753c2b3bed&chksm=9b3c8235ac4b0b238c4902ad7b3944b285baaa0a7ec6f5c7d08a596fdf89735dc7a7afab91f0&scene=178&cur_album_id=2802463431357661184#rd)  
# [文献阅读（十八）、紧跟生信"钱"沿,胰腺癌&免疫多模态图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247489045&idx=1&sn=29feb00f094a9a7e20418b22aa54dd8c&chksm=9b3f6f45ac48e653d4223838eb1d96e65e5ee51ac3f7782841d9b33013bdb67b5e22e9599fac&scene=21#wechat_redirect)  
# [文献阅读（十九）、原发头颈癌和肿瘤转移微生态](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247489511&idx=1&sn=76ec706ebd97e3edd08687f5cf3b6816&chksm=9b3f6eb7ac48e7a10bb81aec5332ecc0544a683f576ea680c8baa621e504816bf65c7d331bdf&scene=21#wechat_redirect)  
# [文献阅读（二十一）、肾脏疾病&代谢&核受体ESRRA](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490575&idx=1&sn=b931a1d50ddcfc4e34999763139695a0&chksm=9b3f655fac48ec49eb2fb01ce8acb15e230091ab3c29ce1471747850be5dae54dbe97fd103b0&scene=21#wechat_redirect)
# [文献阅读（二十二）、PD-1&急性髓细胞性白血病&T Cell](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490929&idx=1&sn=82bf93ccbd9dd906caad8ec709c86153&chksm=9b3f6421ac48ed37117794e7c5b878fa78760627e8bced9c37a2c7c67fd8afd58d0cf0add17e&token=1760636281&lang=zh_CN#rd)  
# [文献阅读（二十三）鼻咽癌微环境&scRNA](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491084&idx=1&sn=e450844d4f102caa88cab6deb630a3f4&chksm=9b3f675cac48ee4a6883559859a21d77256479bb8d030c54a91240827d85b925066bc4f8ebb6&token=333450324&lang=zh_CN#rd)  
# [文献阅读(二十四)、《Nature》：MYB调控衰竭性T细胞对检查点抑制的响应](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491909&idx=1&sn=5d47d2917397248dd6a3a9e85f2acd6d&chksm=9b3c9815ac4b11037e52ba00215a0b6cd23b554139955123ef21ac2cd07a646e97309e7def78&scene=21&token=7173752&lang=zh_CN#wechat_redirect)  
# [文献阅读(二十八)、单细胞都能活检测序了？](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247493604&idx=1&sn=d0f8738e26a4b1c870d329b6a4a7b689&chksm=9b3c9eb4ac4b17a2d0b04d72c350c20f575701de7179f77fcd3567f4b1533cf35df4157f67e1&token=416896563&lang=zh_CN#rd)  
# [文献阅读(二十九)、单细胞测序做到什么程度能毕业|硕士篇](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247493808&idx=1&sn=52c7d330ffdda47c3409d36fe9cbbcd8&chksm=9b3c91e0ac4b18f6edf59d171a0d982b1e3f842148046d8eed19343484feaa84dbd5b7682a11&token=416896563&lang=zh_CN#rd)  
# [文献阅读(三十)、单细胞测序做到什么程度能毕业|博士篇](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494017&idx=1&sn=a2dc4e2eef2bcdb40f5074b665dfa45f&chksm=9b3c90d1ac4b19c79ac5128f96d27789241357c9745bb0389ff5f79f96ab825d3e1c9b29e3af&token=416896563&lang=zh_CN#rd)  
# [2022年了，都有哪些器官/组织有scRNA-Seq数据|小鼠篇](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494796&idx=1&sn=92c91a65719bdf39658814eb691a780f&chksm=9b3c95dcac4b1cca5547440a5881e650da3ea430bad3e6017cbd2891606b2ec40e3ec343bc0e&scene=21&token=468535359&lang=zh_CN#wechat_redirect)  
# [2022年了，都有哪些器官/组织有scRNA-Seq数据|人类篇](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247494806&idx=1&sn=1d15ad9cb4d350a164d3b941bc715233&chksm=9b3c95c6ac4b1cd054cedbd3e8c8f57ee88cc59a6ebda2fd11ba8fea79f5e9e8c65d21ff338c&scene=21&token=468535359&lang=zh_CN#wechat_redirect)  
# [2023年了, 都哪些物种有scRNA-Seq数据](https://mp.weixin.qq.com/s/otvcVXJDMbF_hdQPaPw5hQ)   
# [终于读到一篇用monocle3做拟时序的文章](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247496064&idx=1&sn=0be38306f9cd8c6b819d58cc2e779b87&chksm=9b3c88d0ac4b01c6773b6c4f1a063788509baa27c57e2098b16b8ebe38cc625a87b9d89272c7&scene=21&token=445940587&lang=zh_CN#wechat_redirect)  
# [单细胞转录组+亚细胞空间代谢组=25分文章](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247496126&idx=1&sn=8cb6040a38eb844bb215c0d395196e3a&chksm=9b3c88eeac4b01f8dcb8f4c7a19c5a58841d5fe1623f0a4c585c6dcc281cbbdd2830a1150f6a&scene=21&token=445940587&lang=zh_CN#wechat_redirect)  
# [基于scRNA-Seq&空间转录组的AKI研究](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503664&idx=1&sn=39519d9cfa54e2d5b757f6c249ccbc86&chksm=9b3cb660ac4b3f768c8cfb708d8355c426b2d4ad39538acd5584aa9741ccdf1076b7112492df&token=436478285&lang=zh_CN#rd)  
# [《Nature Methods》: NiCheNet, 能够考虑到胞内转录调控的细胞通讯软件](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503617&idx=1&sn=18389ffe1d201b1c0592ebc923aac851&chksm=9b3cb651ac4b3f4777cb0a6282ad99cd2aad711050d1590f20f970d843f05fb42f8cc02833e2&scene=178&cur_album_id=2802463431357661184#rd)      
# [数百万级单细胞数据时代的多组学分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503603&idx=1&sn=c146eec4528f418226f0103030cb803c&chksm=9b3cb7a3ac4b3eb5bd8be7e41c8c440c80e8223fb10634bfc6d523fe389e3aa0bc14e5ea87d5&scene=178&cur_album_id=2802463431357661184#rd)  
# [两万字长文|当前单细胞 RNA-seq 分析的最佳实践](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503544&idx=1&sn=d563497cc912950d07c2284709f3aede&chksm=9b3cb7e8ac4b3efe66fa4e743f9d6c245a7dd6ea370cd0ee397e92761a2d680cac5c6dd4ea87&scene=178&cur_album_id=2802463431357661184#rd)  
# [《Nature》: 单细胞解析大脑感知流感过程](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503508&idx=1&sn=72781ff97ecae1c32815fe2c4d71b24e&chksm=9b3cb7c4ac4b3ed2b46c8c22f4ce8f93e3f93004ee76542bdd92039d301cec4a10077bb973db&scene=178&cur_album_id=2802463431357661184#rd)  
# [scRNA-Seq+WGCNA+铁死亡：IF=9.69](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503498&idx=1&sn=fcf38bd6fd1dbaf80e203926bb63d5c4&chksm=9b3cb7daac4b3ecca9f2dd590fc581f5b1ae67a7f89be6f76eb438a017b5c335d8a01830b6a5&scene=178&cur_album_id=2802463431357661184#rd)  
# [<IF=27.4> IgG4相关性疾病中颌下腺和血液中的细胞和分子改变](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503081&idx=1&sn=29890957ab75f050d466ca659292b529&chksm=9b3cb5b9ac4b3cafe279783c49cf7928819c5caeffbc34e52ae0c79f5e701ce8a25dd6f0c913&scene=178&cur_album_id=2802463431357661184#rd)    
# [急性心肌梗塞(AMI)外周血scRNA-Seq分析流程及结果梳理](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500943&idx=1&sn=c85cf5be55fa2ae4b85c70ca7e4dfffc&chksm=9b3cbddfac4b34c99abb307f71883ca495e2ed876556a7ad064a3cd69e17827dc534317b372e&scene=178&cur_album_id=2802463431357661184#rd)  
# [酸了，六个样本的scRNA-Seq+差异分析=9分文章](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500273&idx=2&sn=c95ac7ceecce0e6364f49e75058caabe&chksm=9b3cb8a1ac4b31b76b36635ac5e3d3385db6cc8c3fa916541bcb61d9a3ad4417f4bfdaa10274&scene=178&cur_album_id=2802463431357661184#rd)  
# [自测scRNA-Seq+scWGBS=3分三区文章？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247498939&idx=1&sn=1b33d219c795f2de8ff66a5c330a3d24&chksm=9b3c85ebac4b0cfd050a7c0abe11dead991a988aa5b790eac6fd57de2f76a85ef2191d8b9ca6&scene=21#wechat_redirect)  
# [百万级单细胞多组学数据集成](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499336&idx=1&sn=864c535b38ce630487fccc602d893a80&chksm=9b3c8718ac4b0e0e8a0e2560a83546167a82616d99503606e3ac0ef29e90b01cc762e95d64e6&scene=21#wechat_redirect)  
# [空间转录组与单细胞转录组整合分析工具大比拼](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247500243&idx=1&sn=b0e513e85a202ec7ca4041a700a7ae10&chksm=9b3cb883ac4b3195cffd0398d880468b12dd068f9a9fd96b7d0f8aed8961be34293e6ff4a2bb&scene=21#wechat_redirect)  
# [<IF=27.4> IgG4相关性疾病中颌下腺和血液中的细胞和分子改变](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503081&idx=1&sn=29890957ab75f050d466ca659292b529&chksm=9b3cb5b9ac4b3cafe279783c49cf7928819c5caeffbc34e52ae0c79f5e701ce8a25dd6f0c913&scene=21#wechat_redirect)  
# [微生物组研究技术——微生物也能做单细胞](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247488278&idx=1&sn=3662cde0f439a7e7cad7fb0d50082c05&chksm=9b3f6a46ac48e350318b38f1bfe252dd4c315ad058dc4b9aa02fadf2be2291fd0dac50f24cc6&scene=21#wechat_redirect)  
# [scRNA-Seq+WGCNA+铁死亡：IF=9.69](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503498&idx=1&sn=fcf38bd6fd1dbaf80e203926bb63d5c4&chksm=9b3cb7daac4b3ecca9f2dd590fc581f5b1ae67a7f89be6f76eb438a017b5c335d8a01830b6a5&scene=21#wechat_redirect)  
# [《Nature》: 单细胞解析大脑感知流感过程](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503508&idx=1&sn=72781ff97ecae1c32815fe2c4d71b24e&chksm=9b3cb7c4ac4b3ed2b46c8c22f4ce8f93e3f93004ee76542bdd92039d301cec4a10077bb973db&scene=21#wechat_redirect)  
# [两万字长文|当前单细胞 RNA-seq 分析的最佳实践](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503544&idx=1&sn=d563497cc912950d07c2284709f3aede&chksm=9b3cb7e8ac4b3efe66fa4e743f9d6c245a7dd6ea370cd0ee397e92761a2d680cac5c6dd4ea87&scene=21#wechat_redirect)  
# [数百万级单细胞数据时代的多组学分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503603&idx=1&sn=c146eec4528f418226f0103030cb803c&chksm=9b3cb7a3ac4b3eb5bd8be7e41c8c440c80e8223fb10634bfc6d523fe389e3aa0bc14e5ea87d5&scene=21#wechat_redirect)  
# [《Nature Methods》: NiCheNet, 能够考虑到胞内转录调控的细胞通讯软件](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503617&idx=1&sn=18389ffe1d201b1c0592ebc923aac851&chksm=9b3cb651ac4b3f4777cb0a6282ad99cd2aad711050d1590f20f970d843f05fb42f8cc02833e2&scene=21#wechat_redirect)  
# [基于scRNA-Seq&空间转录组的AKI研究](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503664&idx=1&sn=39519d9cfa54e2d5b757f6c249ccbc86&chksm=9b3cb660ac4b3f768c8cfb708d8355c426b2d4ad39538acd5584aa9741ccdf1076b7112492df&scene=21#wechat_redirect)  
# [自测蛋白质组+Bulk RNA-Seq & 挖掘scRNA-Seq=13.6分JASN](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247504408&idx=1&sn=da427a6687b06a82d9dd58ae9f9aa474&chksm=9b3cab48ac4b225e8272d0475ae9295c1ea96172419aa604526f1b8a898748cb3897508d95a0&scene=21#wechat_redirect)  
# [非监督式聚类在scRNA-Seq应用中的挑战](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247505257&idx=1&sn=a524a42b966007defbb2529d300d4c28&chksm=9b3cac39ac4b252f9f88181b932fdc279507c2affc765b53dbccf73084f86f0de25ab1af39bc&scene=21#wechat_redirect)  
# [Seurat V5的《Nature Biotechnology》](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247505826&idx=1&sn=5b9e4aeb41fc2f2943aeac9e22d4d038&chksm=9b3caef2ac4b27e4d57e75a508f00b5605ab99a7b51024110f54778683e14ed8bbcd743f0588&scene=21#wechat_redirect)  
# [当空转遇上单细胞~](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247505721&idx=1&sn=19d02b3858ebcf483d9e071b8ba00092&chksm=9b3cae69ac4b277f5b90fd75819ac8435226ad6daaf744f367d49ceb8099408e097d48656e2f&scene=21#wechat_redirect)  
# [人类肾脏参考组织图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506567&idx=1&sn=e895c0d5a6dde49be76429e9961e61bf&chksm=9b3ca3d7ac4b2ac15dc4e2665e5d624b92bf85d8c94ac70c26e476444451c40a5286c853628f&scene=21#wechat_redirect)  
# [有免疫特征的基质细胞，有没有可能是双细胞？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506493&idx=1&sn=ed65c6b0738a1fb8492da0b4ecdc9c63&chksm=9b3ca36dac4b2a7b4ca45f169526c6ed51777e20811d50c2892e683eee405c25607a42cc7a66&scene=21#wechat_redirect)    
# [CAR-T治疗急性髓细胞性白血病的scRNA-Seq图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506371&idx=1&sn=a9f44bdfd795a4d04615bc307cd23cf5&chksm=9b3ca093ac4b29853f7f068e3969384ed07cd2cf4d64a55c423cd74bce83b914c96844a4cb3d&scene=21#wechat_redirect)  
# [临床样本+单细胞+空转=IF4.8？](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247508744&idx=1&sn=1bc02443acab6d1720fe64e0334ce57a&chksm=9b3cda58ac4b534e3bf7ed2ab225d104d9976c5dc21b8443f88d3a555c5db67449b622ac91b4&scene=21#wechat_redirect)  
# [鼻咽癌的Bulk RNA-Seq与scRNA-Seq联合分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247510548&idx=1&sn=58048c46e837b1d56dffa282272c6d96&chksm=9b3cd344ac4b5a52369e91b4e61fa0a2c951b82cc21165b99f42ba2c2191d4add803bd73c928&scene=21#wechat_redirect)  
# [GPTCelltype：玩GPT4玩出Nature](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247510790&idx=1&sn=b56acf349e417f4ce1a09485256a5f87&chksm=9b3cd256ac4b5b400af5a809c5f595c8e1442513d277501c5bca1c80ce85f18a06a1d8994a9c&scene=21#wechat_redirect)  
# [《Nature Methods》教你如何挑选空转平台](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247520891&idx=1&sn=147d185643bd5dbba4e58b20b2eed872&chksm=9b3ceb2bac4b623d0273eb53b44709e951a493d79338d622b07e256cb4bf385c7d77983e7f48&scene=21#wechat_redirect)   
# [一文了解sc与sn的异同](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247526966&idx=1&sn=9b2aaa1059e6ed68ada65fe004190108&scene=21#wechat_redirect)
#                                          
# ## 单细胞注释复写
# [单细胞注释复写(一):Human Fetal Kidney](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490504&idx=1&sn=0066b2b33067cedb9d1c28689a53faab&chksm=9b3f6298ac48eb8e12cdc7bdd6263e5feb9039f549f7b9dbf4e258854ee61476d09f702f1447&scene=21#wechat_redirect)  
# [单细胞注释复写(二):Human colorectal](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247491642&idx=1&sn=2246fb702277ddc690f056fc8854142e&chksm=9b3c996aac4b107c059a403978eca7be1ae7987659626272c58430910ad75892ad470517aefb&scene=21#wechat_redirect)  
# [单细胞复写|3.急性心肌梗塞(AMI)外周血scRNA-Seq分析实战(链接重置)](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247501145&idx=2&sn=12c2e8d104ad343fc13acb8cda266e6c&chksm=9b3cbc09ac4b351f4ed708a0ed577a21f48ccb0d4173aed94a027b4170147ab85931ee51d1e5&payreadticket=HFvx9c1UB2uxoz-jC4onxEBsjXz1koi2pUYaGWVOLg3rIO1tAGAQrpPmfpXluNpTIBGUd4M#rd)  
# [单细胞复写|4. GSE157783数据集复现及差异分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247503585&idx=1&sn=0035977a80621f09d0774d3c9b0a86cd&chksm=9b3cb7b1ac4b3ea76685ed44abcfb031c698835d9b9b42b269496f6104cd8c9d45467de7eddd&token=436478285&lang=zh_CN#rd)  
# [《NC》代码复现| scRNA-seq揭示人纤维化皮肤病中的成纤维细胞异质性和间充质成纤维细胞增加](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247525649&idx=1&sn=6145f50422ef8680be2a7ffd6c66f5f3&scene=21#wechat_redirect)   
# [代码复写《Developmental Cell》|利用scRNA-seq对新生儿~成人肝脏发育和成熟的时间分析](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247526426&idx=1&sn=7191649a25332d1b8bf2bd9c797d06a7&scene=21#wechat_redirect)   
# [《Nature Communications》| scRNA-seq揭示头颈癌早期转移过程中潜在的免疫逃避机制](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247526426&idx=2&sn=e33c59e785800a498b57d7a971a6e27a&scene=21#wechat_redirect)   
# [《Cell》复现|中性粒细胞scRNA-seq图谱阐明抗肿瘤抗原呈递效力](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247526776&idx=1&sn=99eb6a16c7bf345445b6ff4b61110138&scene=21#wechat_redirect)  
# [复现24年前的《Nature》—GEO中的第一个数据](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247526262&idx=1&sn=0cf9dc3dfb55a636634612a2aa1d2e23&scene=21#wechat_redirect)  
# [CellRank的教程重现](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507012&idx=1&sn=fb01cca63326fd6e8c682434a53fc2ea&chksm=9b3ca514ac4b2c02115ae5dea368190ce74997ac16a07fdc29f76cef815600b1e5b813c71163&scene=21#wechat_redirect)  
# [人类肾脏scRNA-seq图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507501&idx=2&sn=dec30ee34e483c3695613b3e41488a78&chksm=9b3ca77dac4b2e6b1734f1fd3d0134468b2ca8bc9f9118edc6e114a37f7ceea7dc4fa41b4556&scene=21#wechat_redirect)  
# [小鼠纹状体单细胞基因变化图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506900&idx=1&sn=83774f3062e42ca0eee1df2cf19144f0&chksm=9b3ca284ac4b2b9221736162ac44280d6184c9f8b0cc13bbcc4c914c14b707fe13a850c6418e&scene=21#wechat_redirect)  
# [T Cell单细胞参考图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506779&idx=1&sn=dc97a74cfac6caeabaa838160385c3b2&chksm=9b3ca20bac4b2b1da51998764b22ac3c08415298a47ded14283337386f18f9548749de70725e&scene=21#wechat_redirect)  
# [MLP模型鉴定免疫细胞发育和瘤内T细胞耗竭的单细胞染色质图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506754&idx=1&sn=3e3e9b4857bb0371b7c1113df792df2e&chksm=9b3ca212ac4b2b043baf0c72684d186a1c065ff7b3d865c308aff4f3483ab5923638fda8d1ff&scene=21#wechat_redirect)  
# [干细胞样CD4+ T细胞与自身免疫性血管炎](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247506684&idx=1&sn=c27703a05d4fffae6a5e9e54b4fa7dfe&chksm=9b3ca3acac4b2abaa1e5d4f4a7e4d0730f2e5c4d72c7b3d3d29a653c3bae6590d90e31908dde&scene=21#wechat_redirect)  
# [小鼠与人类的空间表观单细胞图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507501&idx=1&sn=94c2a38b98957a9526d510798f791419&chksm=9b3ca77dac4b2e6b40d37e73f54791d55a919f50992bf62ac557ecbac919fb6415dca468dba4&scene=21#wechat_redirect)  
# [人类肾脏scRNA-seq图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507501&idx=2&sn=dec30ee34e483c3695613b3e41488a78&chksm=9b3ca77dac4b2e6b1734f1fd3d0134468b2ca8bc9f9118edc6e114a37f7ceea7dc4fa41b4556&scene=21#wechat_redirect)  
# [胰腺导管腺癌的内异质性和恶性进展scRNA-seq图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507450&idx=1&sn=ee9b7d482ad9a88c6da18341c7ec5254&chksm=9b3ca4aaac4b2dbc83fc6cc4c084636fb397e56d0057dfb5b546acf402661ab64d241f985073&scene=21#wechat_redirect)  
# [肝脏发育、成熟的单细胞图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507353&idx=1&sn=f2922ce2337c31a781a690ea59c5abde&chksm=9b3ca4c9ac4b2ddfc267838fe6d794027643e89fe054a5ca4d48e51976ea825af7fed285e26b&scene=21#wechat_redirect)  
# [微环境中肿瘤免疫屏障决定了免疫治疗的效果](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507337&idx=1&sn=e23d849663a6cfcad19bacdaeeb85b5f&chksm=9b3ca4d9ac4b2dcf8ce4d9720bf7588a30e932adea4297862927ebfe20c29f8961436f761b8c&scene=21#wechat_redirect)  
# [scRNA揭示T1D骨髓与骨质减少的关系](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507144&idx=1&sn=ed6e95f16601e6196cea7179b88db250&chksm=9b3ca598ac4b2c8ecfd1025246eb79c6aa8bf5703f8698d23f1d1f673121e3bf8168a32dba07&scene=21#wechat_redirect)  
# [肝癌单细胞分群marker](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507136&idx=2&sn=16f8afe29fb59b9a8e9cac265098709c&chksm=9b3ca590ac4b2c86c9554bf20edd587f5d758636f55cf9e390220249d3338b1fde7717fb3c6e&scene=21#wechat_redirect)  
# [胶质瘤相关巨噬细胞scRNA图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507096&idx=1&sn=98e42a5768074d2766f2ce64214ac391&chksm=9b3ca5c8ac4b2cdebe68a8519a087e41a4c88a7e465c5ddc107f633a5fbf3fcc8bbff6aa8821&scene=21#wechat_redirect)  
# [人类乳腺癌scRNA-seq数据](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507792&idx=2&sn=52dc1c2e5a23910ca0d8330ea917b60d&chksm=9b3ca600ac4b2f16a04624c74f2bcd3d24d300a835053a6f76ac4e1421b90c34b2ea01f1e7c8&scene=21#wechat_redirect)  
# [肢端黑色素瘤细胞scRNA-seq图谱&免疫治疗](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507783&idx=1&sn=6ca5dc13a6aefbe78cce3f938fe04e40&chksm=9b3ca617ac4b2f01fb6e4677668f052f6361282ff394745394d363b54ea4c1eb58d2f30da30a&scene=21#wechat_redirect)  
# [皮质神经注释复写](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507783&idx=2&sn=d1e28f5f4bddb3a3e4f1a85db1feac7f&chksm=9b3ca617ac4b2f01b40ae6adbe48e9371e50f9df543638e8c3e8e57bf9312b5daf731041ee16&scene=21#wechat_redirect)  
# [非小细胞肺癌免疫治疗后肿瘤微环境重塑scRNA图谱](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507757&idx=1&sn=1045b58a64b4602881e9853d2a908012&chksm=9b3ca67dac4b2f6b4d33f10e71da2a6287968512677e77017b3820bf9aedfcd60301bb14d401&scene=21#wechat_redirect)  
# [小鼠脑动脉瘤模型的scRNA-Seq分析](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507717&idx=2&sn=1d99c171cbfedc11f11fc827b643c89d&chksm=9b3ca655ac4b2f4360790bde80bc1452e9e314c6a0ae2dc6f9d19edbcc580e81ace3080dbb5f&scene=21#wechat_redirect)  
# [儿童白血病的单细胞数据读取](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507693&idx=2&sn=8d45c680ad78f52043d4ce4cbbf3c67a&chksm=9b3ca7bdac4b2eab2405773b5d71ff1f5495250ef813b24429e1b920c1b96a4ee391507d011b&scene=21#wechat_redirect)  
# [疟疾中免疫细胞亚群的调控响应](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247507659&idx=2&sn=a6860bb0eb8d3854a37acd0342777956&chksm=9b3ca79bac4b2e8dc0c3aa83c000894a53e2b316388495703e58d42ab205dd093bf578807b74&scene=21#wechat_redirect)  
# [单细胞空间转录组在脑科学领域的应用](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247520599&idx=1&sn=f5de77dfcb186b29ff0f0657ea680baf&chksm=9b3ce807ac4b61119cff06d6dc1517e05e6296b2039d8a150121964eb134d30259e5f9e965ff&scene=21#wechat_redirect)  
# [snATAC联合scRNA揭示小鼠肾脏细胞分化程序和疾病靶点](http://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247522324&idx=1&sn=62993c2354af1e793fcfd0e06586af17&chksm=9b3ce144ac4b6852aa1c38a44041ccf1bbf3353cfe515bcb28b434a184368e5888641690fd5b&scene=21#wechat_redirect)  

# %% [markdown]
# # 欢迎致谢
# 如果以上内容对你有帮助，欢迎在文章的Acknowledgement中加上这一段，联系客服微信可以发放奖励：
#
# Since Biomamba and his wechat public account team produce bioinformatics tutorials and share code with annotation, we thank Biomamba for their guidance in bioinformatics and data analysis for the current study. 
