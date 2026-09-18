"""infercnvpy 独立 CNV 入口：注释后 h5ad + 正常 reference → 相对 CNV。"""

from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import infercnvpy as cnv
import pandas as pd
import scanpy as sc


# 功能：把当前 Matplotlib 图同时保存为 PDF 和 300 dpi PNG。
def save_current_figure(stem: Path, width: float = 10, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


# 功能：读取主流程表达对象并选择用于 infercnvpy 的表达层。
def read_expression(
    input_h5ad: str,
    layer: str | None,
    use_raw: bool,
    normalize: bool,
):
    adata = sc.read_h5ad(input_h5ad)
    if use_raw:
        if adata.raw is None:
            raise ValueError("--use-raw 已启用，但 h5ad 中没有 adata.raw。")
        raw = adata.raw.to_adata()
        raw.obs = adata.obs.copy()
        adata = raw
    elif layer:
        adata.X = adata.layers[layer].copy()
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


# 功能：从本地 GTF 或 BioMart 添加染色体坐标；本地 GTF 更适合固定基因组版本。
def add_genomic_positions(
    adata,
    gtf: str | None,
    gtf_gene_id: str,
    gene_id_col: str,
    species: str,
    biomart_gene_id: str,
) -> str:
    adata.var[gene_id_col] = adata.var_names
    if gtf:
        cnv.io.genomic_position_from_gtf(
            gtf,
            adata,
            gtf_gene_id=gtf_gene_id,
            adata_gene_id=gene_id_col,
        )
        return str(Path(gtf).expanduser().resolve())
    cnv.io.genomic_position_from_biomart(
        adata=adata,
        adata_gene_id=gene_id_col,
        biomart_gene_id=biomart_gene_id,
        species=species,
    )
    return f"BioMart: {species}/{biomart_gene_id}"


# 功能：运行 infercnvpy、CNV 聚类和 score，并保存审计表及完整对象。
def run_infercnvpy(
    input_h5ad: str,
    output_dir: str,
    object_out: str | None,
    celltype_col: str,
    references: list[str],
    layer: str | None,
    use_raw: bool,
    normalize: bool,
    gtf: str | None,
    gtf_gene_id: str,
    species: str,
    gene_id_col: str,
    biomart_gene_id: str,
    window_size: int,
    leiden_resolution: float,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = (
        Path(object_out).expanduser().resolve()
        if object_out
        else output / "infercnvpy.h5ad"
    )
    output.mkdir(parents=True, exist_ok=True)
    object_path.parent.mkdir(parents=True, exist_ok=True)

    adata = read_expression(input_h5ad, layer, use_raw, normalize)
    if celltype_col not in adata.obs:
        raise KeyError(f"obs 中缺少 reference 列：{celltype_col}")
    adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")
    missing = sorted(set(references).difference(adata.obs[celltype_col].cat.categories))
    if missing:
        raise ValueError(f"Reference categories not found: {', '.join(missing)}")

    coordinate_source = add_genomic_positions(
        adata,
        gtf,
        gtf_gene_id,
        gene_id_col,
        species,
        biomart_gene_id,
    )
    position_columns = [
        column for column in ["chromosome", "start", "end"]
        if column in adata.var
    ]
    adata.var[position_columns].reset_index(names="gene").to_csv(
        output / "genomic_positions_used.csv",
        index=False,
    )

    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color=celltype_col, show=False)
        plt.gca().set_aspect("equal")
        save_current_figure(output / "input_umap_celltype", width=8, height=7)

    # reference_cat 是已知正常细胞；窗口平滑得到的是相对于 reference 的表达偏移。
    cnv.tl.infercnv(
        adata,
        reference_key=celltype_col,
        reference_cat=references,
        window_size=window_size,
    )
    cnv.pl.chromosome_heatmap(adata, groupby=celltype_col, show=False)
    save_current_figure(output / "cnv_chromosome_heatmap_celltype", width=12, height=8)

    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata, resolution=leiden_resolution)
    cnv.pl.chromosome_heatmap(
        adata,
        groupby="cnv_leiden",
        dendrogram=True,
        show=False,
    )
    save_current_figure(output / "cnv_chromosome_heatmap_leiden", width=12, height=8)
    cnv.tl.umap(adata)
    cnv.tl.cnv_score(adata)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    cnv.pl.umap(adata, color="cnv_score", ax=axes[0], show=False)
    cnv.pl.umap(adata, color=celltype_col, ax=axes[1], show=False)
    for axis in axes:
        axis.set_aspect("equal")
    save_current_figure(output / "cnv_score_umap", width=12, height=5.5)

    columns = [celltype_col, "cnv_leiden", "cnv_score"]
    metadata = adata.obs[columns].copy()
    metadata.insert(0, "cell", metadata.index)
    metadata.to_csv(output / "infercnvpy_metadata.csv", index=False)
    adata.write_h5ad(object_path)

    (output / "reference_categories.txt").write_text(
        "\n".join(references) + "\n",
        encoding="utf-8",
    )
    (output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# infercnvpy 结果说明",
                "",
                f"- 输入 h5ad：{Path(input_h5ad).expanduser().resolve()}",
                f"- 表达来源：layer={layer}；use_raw={use_raw}；normalize={normalize}",
                f"- reference：{', '.join(references)}",
                f"- 基因坐标：{coordinate_source}",
                f"- window size：{window_size}；CNV Leiden resolution：{leiden_resolution}",
                f"- infercnvpy 对象：{object_path}",
                "- reference 身份和基因组坐标版本是关键分析假设，必须与物种和 gene ID 匹配。",
                "- CNV score/cluster 是转录组相对 CNV 推断，不能替代 DNA 拷贝数检测。",
                "- 不对单细胞 CNV score 做 CON/EXP 显著性检验；条件比较应使用样本级汇总。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (output / "environment_versions.txt").write_text(
        "\n".join(
            [
                f"scanpy=={version('scanpy')}",
                f"infercnvpy=={version('infercnvpy')}",
                f"anndata=={version('anndata')}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


# 功能：定义 h5ad、reference、表达层和基因坐标参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run infercnvpy from an annotated h5ad.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--object-out", default=None)
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument("--references", required=True, help="Comma-separated normal cell types.")
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use-raw", action="store_true")
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--gtf", default=None)
    parser.add_argument("--gtf-gene-id", choices=["gene_id", "gene_name"], default="gene_name")
    parser.add_argument("--species", default="hsapiens")
    parser.add_argument("--gene-id-col", default="gene_ids")
    parser.add_argument("--biomart-gene-id", default="external_gene_name")
    parser.add_argument("--window-size", type=int, default=250)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_infercnvpy(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        object_out=args.object_out,
        celltype_col=args.celltype_col,
        references=[item.strip() for item in args.references.split(",") if item.strip()],
        layer=args.layer,
        use_raw=args.use_raw,
        normalize=args.normalize,
        gtf=args.gtf,
        gtf_gene_id=args.gtf_gene_id,
        species=args.species,
        gene_id_col=args.gene_id_col,
        biomart_gene_id=args.biomart_gene_id,
        window_size=args.window_size,
        leiden_resolution=args.leiden_resolution,
    )