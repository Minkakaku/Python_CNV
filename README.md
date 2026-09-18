# R-Python-infercnv-copykat-workflow-template

独立 RNA-CNV 分析仓库。每种方法对应一个根目录完整脚本，直接接收主流程输出：

- infercnv_workflow.R：Seurat RDS + gene-order 文件，运行 Broad inferCNV。
- copykat_workflow.R：Seurat RDS，按样本分别运行 CopyKAT。
- infercnvpy_workflow.py：AnnData h5ad + GTF/BioMart 坐标，运行 infercnvpy。

## inferCNV

    Rscript infercnv_workflow.R --input annotated.rds --gene-order gene_order.tsv --output results/infercnv --references "T cell,B cell,Endothelial cell"

输入必须使用 raw counts。gene-order 文件必须与物种、gene ID 和基因组版本一致；正常 reference 需要由 marker 和样本来源支持。

## CopyKAT

    Rscript copykat_workflow.R --input annotated.rds --output results/copykat --sample-col orig.ident --references "T cell,B cell"

CopyKAT 官方建议每个样本单独运行，合并样本可能把批次差异误判为 CNV。脚本会按 sample-col 自动拆分。human 使用 genome=hg20，mouse 使用 genome=mm10；小鼠预测需要额外目视检查。

## infercnvpy

    python infercnvpy_workflow.py --input-h5ad annotated.h5ad --output-dir results/infercnvpy --celltype-col cell_type --references "T cell,B cell" --gtf genes.gtf

本地 GTF 更便于固定基因组版本。若输入层是 raw counts，可增加 --layer counts --normalize；若 X 已为 log-normalized expression，则不增加 --normalize。

## 解释边界

三种方法都根据 RNA 表达推断相对 CNV，不能替代 DNA 测序。正常 reference、基因组坐标和物种是关键分析假设。不要把单细胞 CNV score 当作生物学重复进行条件显著性检验。

legacy/original/ 保存旧项目脚本，source_catalog/ 保存来源语义副本，均不作为正式入口。

方法依据：inferCNV 官方 Wiki、infercnvpy 官方文档和 CopyKAT 论文（PMID: 33462507）。