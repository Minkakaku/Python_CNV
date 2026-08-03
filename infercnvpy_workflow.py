"""Portable infercnvpy workflow converted from the specific inferCNV notebook."""

from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import infercnvpy as cnv
import scanpy as sc


def save_current_figure(stem: Path, width: float = 10, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_infercnvpy(
    input_h5ad: str,
    output_dir: str,
    object_out: str,
    celltype_col: str,
    references: list[str],
    species: str,
    gene_id_col: str,
    biomart_gene_id: str,
    window_size: int,
    leiden_resolution: float,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = Path(object_out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    object_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_h5ad)
    adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")
    missing_references = sorted(set(references).difference(adata.obs[celltype_col].cat.categories))
    if missing_references:
        raise ValueError(f"Reference categories not found: {', '.join(missing_references)}")

    adata.var[gene_id_col] = adata.var_names
    cnv.io.genomic_position_from_biomart(
        adata=adata,
        adata_gene_id=gene_id_col,
        biomart_gene_id=biomart_gene_id,
        species=species,
    )

    sc.pl.umap(adata, color=celltype_col, show=False)
    plt.gca().set_aspect("equal")
    save_current_figure(output / "input_umap_celltype", width=8, height=7)

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
    sc.pl.umap(adata, color="cnv_score", ax=axes[0], show=False)
    sc.pl.umap(adata, color=celltype_col, ax=axes[1], show=False)
    for axis in axes:
        axis.set_aspect("equal")
    figure.tight_layout()
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
                "# infercnvpy 结果",
                "",
                f"- 输入对象：{Path(input_h5ad).expanduser().resolve()}",
                f"- 细胞类型列：{celltype_col}",
                f"- reference：{', '.join(references)}",
                f"- 物种：{species}",
                f"- window size：{window_size}",
                f"- CNV Leiden resolution：{leiden_resolution}",
                f"- 结果目录：{output}",
                f"- infercnvpy 对象：{object_path}",
                "- CNV score 和 cluster 是转录组推断结果，不能替代 DNA 层面的拷贝数检测。",
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run infercnvpy from an annotated h5ad file.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument(
        "--references",
        required=True,
        help="Comma-separated normal/reference cell-type labels.",
    )
    parser.add_argument("--species", default="mmusculus")
    parser.add_argument("--gene-id-col", default="gene_ids")
    parser.add_argument("--biomart-gene-id", default="external_gene_name")
    parser.add_argument("--window-size", type=int, default=250)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument(
        "--output-dir",
        default="results/06_infercnv_copykat/infercnvpy",
    )
    parser.add_argument(
        "--object-out",
        default="data/processed/infercnvpy.h5ad",
    )
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_infercnvpy(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        object_out=args.object_out,
        celltype_col=args.celltype_col,
        references=[item.strip() for item in args.references.split(",") if item.strip()],
        species=args.species,
        gene_id_col=args.gene_id_col,
        biomart_gene_id=args.biomart_gene_id,
        window_size=args.window_size,
        leiden_resolution=args.leiden_resolution,
    )
