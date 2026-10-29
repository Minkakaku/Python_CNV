# GSE160269 平台原码 infercnvpy 结果

本结果完成于 2026-10-28，2026-10-29 归档至此。51,571 个细胞来自 64 个样本，其中上皮 43,834、巨噬细胞参考 7,737。使用全队列共享参考池合并运行一次，不是逐样本结果。

## 分析与注释

输入为 27,431 个基因的原始 counts，经每细胞 10,000 总量归一化及 log1p 后用于 infercnvpy。公司平台原始计算、人工注释和 Rmd 模板文件未修改，来源哈希见上一级说明。

infercnvpy 0.4.3、scanpy 1.9.5；GRCh38；`reference_key='cell_type'`、`reference_cat=['Macrophage']`、`window_size=250`；保留原默认 PCA 50、neighbors 15、Leiden resolution 1、UMAP min_dist 0.5。得到 2,354 个 CNV 窗口和 60 个 CNV 簇。输入按固定种子 20261028 打乱行顺序，使默认 5,000 细胞技术块混合上皮与参考，未按样本拆分模型。

RNA UMAP 在共同的 Harmony 30 维坐标上重新计算，包含上皮和巨噬细胞；保留 Epi1–Epi16。CNV 的 60 个簇与上皮转录亚群是不同的分组。

原平台需要人工注释。用户批准的保守映射：normal 为 2、3、6、13、20、26、30、53；uncertain 为 1、7、8、34、56；其余为 tumor 候选。全体细胞 normal 10,333、uncertain 6,637、tumor 34,601。得分为相对 CNA 强度，不是恶性概率或整数 DNA 拷贝数。

## 文件

- [report.html](report.html)：原模板生成的原始 HTML，内含七张科学图和完整可检索细胞表。
- [CNV/](CNV/)：11 对原生 PNG/PDF；`metadata.csv` 是交付表，状态列名为 `cnv_status`；`meta_inferCNV.xls` 是原报告使用的制表符文本，状态列名为 `inferCNV`，并非 Excel 二进制格式。两列的细胞状态完全相同。
- [supplement/](supplement/)：Ref/Epi1–Epi16 评分图及实际入图数据、逐簇审阅表、用户批准映射、正常组织对照、原核验记录及计算环境版本。

报告、图表和原有表格均原样复制。此前的错误版本、CopyKAT 和 64 个独立样本分析不属于本目录。本次没有重新运行分析。

参考细胞的细胞类型表达差异、样本和批次效应会影响 RNA-CNV。tumor 是人工映射下的候选；技术核验不等同于独立生物学验证，本结果未经 DNA CNV 或病理真值验证。
