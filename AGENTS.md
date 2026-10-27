# 项目说明

- 根目录：`/home/lcpyspa/spa`。除 Codex 配置与 Conda 环境操作外，所有文件操作限于本项目。
- 远端：`https://github.com/LouChen-med/spatial_learn.git`。
- 内容：空间组学学习 Notebook，按学习单元保存在 `spa04/` 与 `spa05/`。
- 环境管理优先使用 `/home/lcpyspa/miniconda3/bin/conda`；本次未修改环境，依赖尚未验证。
- 保留现有 Notebook 及输出；所有层级的 `.ipynb_checkpoints/` 不纳入 Git。
- 2026-10-27：新增 `spa04/ST_introduction_lc.py`，按原顺序转换全部 125 个单元（50 个代码、75 个 Markdown），代码保持不变，非代码内容转为注释并用 `#%%` 分隔；逐单元核对内容，未运行分析，未修改 Conda 环境或 Git。
- 每个子目录维护独立 `README.md`；新增代码、输入、输出、脚本和临时文件应按用途分目录。
- 代码保持精炼、清晰、易维护。Conda 或 Git 修改命令记录在项目内 `.txt` 文件，包含工作目录。
- 2026-10-27：已初始化 main 分支并提交两个 Notebook 与文档；推送因缺少 GitHub HTTPS 凭据失败，尚未上传。命令见 `git_operations.txt`，配置认证后重试 `git push -u origin main`。
- 2026-10-27：修复 `spa04/ST_introduction_lc.py` 的 Notebook 导出语法：将 Bash 单元及 `!` 命令改为 `subprocess` 调用，裸日志和说明文字转为注释，修复 Shell 续行符后的空格；保留单元顺序、教学内容及原始 Notebook，示例路径未改为本机路径。仅做静态语法检查，不执行安装、下载或分析，不修改 Conda 环境或 Git。
- 2026-10-27：将 `spaceranger count` 及结果查看单元的输出固定到 `/home/lcpyspa/spa/spa04/Visium_HD_Mouse_Brain`，使用 `--output-dir`，输入通过 `DATA_DIR` 统一引用；原 `/data/04_ST_introduction/data/spaceranger` 为用户提供但尚未核实的路径，未越出项目访问。检查时项目内无分析输入，目标目录为空；移除后台 `&`，未运行分析、下载数据或修改 Conda/Git。原始 Notebook 保持不变。
- 2026-10-27：按用户要求提交当前 Notebook、Python 脚本及文档，保留用户已有修改；确认 `.gitignore` 排除全部 Notebook 检查点、Space Ranger 运行目录、FASTQ 目录及约 27 GB 的 FASTQ 压缩包。两个 Notebook 的 JSON、Python 脚本语法及待上传文件的常见凭据模式检查通过；未运行分析或修改 Conda。远端 `origin/main` 读取预检成功，但推送认证预检因缺少凭据失败，实际推送遇到 GnuTLS 连接中断，HTTP/1.1 重试 60 秒后超时；本地提交保留，尚未完成上传。本次 Git 修改命令及执行结果记录在 `git_operations.txt`，需恢复连接并配置 GitHub 认证后重试推送。
