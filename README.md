# R-Python-infercnv-copykat-workflow-template

inferCNV 与 CopyKAT 肿瘤拷贝数分析代码素材模板。

## 正式 Python 入口

`specific/infercnv.ipynb` 已转换为根目录的 `infercnvpy_workflow.py`，正式分析不再使用 notebook。reference cell types、物种、窗口和输出路径均通过参数指定：

```bash
python infercnvpy_workflow.py \
  --input-h5ad data/processed/annotated.h5ad \
  --celltype-col cell_type \
  --references Macrophage,T,DC,B,Endothelial \
  --species mmusculus
```

默认输出：

- 图表和 metadata：`results/06_infercnv_copykat/infercnvpy/`
- 完整对象：`data/processed/infercnvpy.h5ad`
- 每张正式图均保存 PDF 和 300-dpi PNG。
- `README_results_CN.md` 会记录 reference、物种、window size 和主要限制。

基因位置默认通过 Biomart 查询；`mmusculus + external_gene_name` 对应原 notebook 的小鼠 gene symbol 用法。人类数据应改为匹配的 species 和 gene identifier。
## 主要代码

- `infercnv/infercnv_seurat.R`：基于 Seurat 对象的 inferCNV 流程。
- `infercnv/examples/`：项目示例。
- `infercnv/report.rmd`：结果报告模板。
- `copykat/copykat_20260317.R`：CopyKAT 分析示例。
- `infercnv/pip_infercnv.py`、`platform_wrappers/`：原平台封装。

## 使用前注意

`infercnv_seurat.R` 会查找单细胞模板中的 `00_functions.R`；在新项目中应从 `R-scrna-workflow-template` 选择相应公共函数。CopyKAT 示例含固定输入对象、reference cell type 和项目参数，运行前必须改为当前项目设置。

原始精选脚本保存在 `../R-Python-bioinformatics-code-archive/legacy/original/06_infercnv_copykat/`；specific notebook 原稿保存在 `../R-Python-bioinformatics-code-archive/legacy/original/specific/infercnv.ipynb`。

## 其他人员来源代码

`source_catalog/` 收录 78 份 inferCNV/CopyKAT 语义副本，并在 `SOURCE_CATALOG.tsv` 中区分 inferCNV 与 CopyKAT、owner、场景和哈希。根目录现有脚本仍作为首选素材，来源版本用于补充不同 reference、物种和平台实现。
