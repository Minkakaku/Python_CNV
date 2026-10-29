# scVision5 inferCNV 报告与结果

按用户要求，于 2026-10-29 将公司原始示例和 GSE160269 的平台原码结果作为结果档案收录。此目录是仓库通常不提交项目产物规则的一次明确例外；实际分析项目仍在本地维护。

| 结果 | 细胞数 | 内容 |
|---|---:|---|
| [公司原始示例](company_demo/) | 49,058 | 用户提供的 `6.inferCNV分析.zip` 的完整解压内容 |
| [GSE160269](GSE160269/) | 51,571 | 全队列合并计算，43,834 上皮 + 7,737 巨噬参考；原报告、全部原生图表及补充审阅材料 |

两套结果各自保存 `report.html` 和 `CNV/`，没有混合细胞、得分、注释或图表。HTML 可以下载后在浏览器打开；GitHub 文件页显示的是 HTML 源码。报告、PNG/PDF 和原有表格均保持来源文件的原始字节，未重新绘图或渲染。

原版 Rmd 的使用入口见 [templates/scvision5](../../templates/scvision5/)。模板、计算和人工注释代码的 SHA-256：

| 来源档案 | SHA-256 |
|---|---|
| [report.rmd](../../legacy/original/infercnv/report.rmd) | `60fbdbcac5337961d328c82185fd9192eb8e9bd4fa4e0a85295fe46a6ab091d1` |
| [pip_infercnv.py](../../legacy/original/infercnv/pip_infercnv.py) | `1c5f8835dab2572d5f4901efa3079d19ef19f9ce6519fe6d8ed54e958da200d1` |
| [run_infercnv.py](../../legacy/original/platform_wrappers/run_infercnv.py) | `0a2de74072ffef011ee2b4770a3b23db5f0f0141da7962faf111a30ae26792cb` |

`SHA256SUMS.tsv` 记录此目录所有结果文件及模板副本的字节数和 SHA-256。公司 ZIP 原件 SHA-256 为 `f5a5092577482296b423bb97388cc6392c62a066e7432973ce40b84af01cd1e1`；为避免重复存储，仓库收录其完整解压文件，不再重复放 ZIP。

GSE160269 是用公司平台原代码运行的公开数据分析，不是公司出具或认可的分析。公司名称保留自原始报告模板。CNV 状态来自 RNA 推断和人工映射，不能替代 DNA CNV 检测。

本次仅归档已有报告、图表、注释与审阅材料，不包含原始 RDS、H5AD、表达矩阵或已停止的 CopyKAT 分析；没有重跑分析。
