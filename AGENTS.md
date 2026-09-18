# RNA-CNV 项目约定

- 本仓库是独立分析点；正式入口仅为根目录三个完整 R/Python 脚本。
- 输入路径、reference、物种、基因坐标和输出目录通过命令行传入。
- inferCNV/CopyKAT 使用 raw counts；infercnvpy 使用非中心化表达并显式记录是否归一化。
- CopyKAT 按生物学样本分别运行。
- 图成对保存 PNG/PDF，表格和中文说明放同一输出目录。
- 所有结论均标明 RNA 推断不能替代 DNA CNV。
- legacy/original/ 与 source_catalog/ 保持原字节。