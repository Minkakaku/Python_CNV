# scVision5 原版 inferCNV Rmd

[report.rmd](report.rmd) 与 [来源档案](../../legacy/original/infercnv/report.rmd) 逐字节相同，SHA-256 为 `60fbdbcac5337961d328c82185fd9192eb8e9bd4fa4e0a85295fe46a6ab091d1`。保留公司署名、正文、CSS、图引用和表格格式。

两套现成报告见 [examples/scvision5](../../examples/scvision5/)，无需重渲染即可查看。

## 对接自己的分析结果

模板读取同目录内以下 PNG（不更改文件名）：

```text
chromosome_heatmap.png
cnv_celltype_bycell.png
cnv_celltype_bycluster.png
cnv_cnvleiden_bycell.png
cnv_cnvleiden_bycluster.png
cnv_vlnplot.png
umap_cnv.png
meta_inferCNV.xls
```

`meta_inferCNV.xls` 是 TSV 文本，第一列为细胞 ID，随后为 `cell_type, cnv_leiden, cnv_score_cell, cnv_score_cluster, inferCNV`。先完成并审核簇注释，再生成此表。不能将平台编辑器初始填充的 normal 当成分析结论。

将 Rmd 和上述输入复制到新的报告工作目录，避免覆盖已交付的报告，然后在具备 R、Pandoc、rmarkdown、knitr、DT 的环境运行：

```bash
Rscript -e "rmarkdown::render('report.rmd', output_file='report.html')"
```

公司原始示例仅提供 `metadata.csv`，未附原生 TSV；不要直接假定它是模板要求的 `.xls` 输入。GSE160269 同时保留交付 CSV 和原生 TSV。

原计算脚本及 Streamlit 入口保存在 `legacy/original/`，含平台固定路径和清理逻辑，是来源档案；本次没有将其改成通用命令行工作流，也没有替换现有 `src/python/infercnvpy_workflow.py`。该原报告模板只适用于其规定的产物和表格，不表示其他工作流的输出可直接套用。
