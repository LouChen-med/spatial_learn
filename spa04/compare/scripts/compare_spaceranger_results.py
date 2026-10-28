#!/usr/bin/env python3

import argparse
import csv
import hashlib
import html
import json
import os
import re
from pathlib import Path

import h5py
import numpy as np
import pyarrow.parquet as parquet


SCRIPT_PATH = Path(__file__).resolve()
SPA04_ROOT = SCRIPT_PATH.parents[2]
PROJECT_ROOT = SPA04_ROOT.parent
RUNS = {
    "local_run": SPA04_ROOT / "Visium_HD_Mouse_Brain",
    "official_download": SPA04_ROOT / "out_10x",
}
SIZES = (2, 8, 16)
PERCENT_METRICS = {
    "Valid Barcodes", "Valid UMI Sequences", "Sequencing Saturation",
    "Q30 Bases in Barcode", "Q30 Bases in Probe Read", "Q30 Bases in UMI",
    "Reads Mapped to Probe Set", "Reads Mapped Confidently to Probe Set",
    "Fraction Reads in Squares Under Tissue",
    "Reads Mapped Confidently to the Filtered Probe Set",
    "Reads Half-Mapped to Probe Set", "Reads Split-Mapped to Probe Set",
    "Estimated UMIs from Genomic DNA", "Fraction of Squares Under Tissue 2 µm",
    "Fraction of Bins Under Tissue 8 µm", "Fraction of Bins Under Tissue 16 µm",
    "Reads in Cells", "UMIs in Cells",
}


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_metrics(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        row = next(csv.DictReader(stream))
    return {"Valid UMI Sequences" if key == "Valid UMIs" else key: value for key, value in row.items()}


def read_summary(path):
    text = path.read_text(encoding="utf-8")
    marker = "const data = "
    offset = text.rfind(marker)
    if offset < 0:
        raise ValueError(f"Cannot find embedded Space Ranger summary data in {path}")
    data, _ = json.JSONDecoder().raw_decode(text[offset + len(marker):])
    commands = []

    def visit(value):
        if isinstance(value, dict):
            if value.get("title") == "Command Line Arguments":
                for row in value.get("data", value.get("card", {}).get("data", [])):
                    for cell in row:
                        for fragment in cell if isinstance(cell, list) else [cell]:
                            commands.extend(re.findall(r"<pre[^>]*>(.*?)</pre>", fragment, re.DOTALL))
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit({key: value for key, value in data.items() if key != "_resources"})
    if not commands:
        raise ValueError(f"Cannot find Space Ranger command line in {path}")
    command = html.unescape(re.sub(r"<[^>]+>", "", commands[0])).strip()
    return {
        "sample": data.get("sample", {}),
        "command": command,
        "alarms": data.get("alarms", {}).get("alarms", []),
    }


def as_text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_features(matrix):
    features = matrix["features"]
    return tuple(as_text(value) for value in features["id"][:]), tuple(
        as_text(value) for value in features["name"][:]
    )


def barcode_grid_indices(barcodes, side):
    chars = barcodes.view(np.uint8).reshape(len(barcodes), barcodes.dtype.itemsize)
    if chars.shape[1] < 21:
        raise ValueError("Expected fixed-width Visium HD square barcodes")
    weights = np.array([10000, 1000, 100, 10, 1], dtype=np.int32)
    rows = (chars[:, 8:13].astype(np.int32) - 48) @ weights
    columns = (chars[:, 14:19].astype(np.int32) - 48) @ weights
    if np.any(rows < 0) or np.any(columns < 0) or np.any(rows >= side) or np.any(columns >= side):
        raise ValueError("Square barcode falls outside the slide grid")
    return rows * side + columns


def matrix_stats(path, side):
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        feature_count, barcode_count = map(int, matrix["shape"][:])
        feature_ids, feature_names = read_features(matrix)
        pointers = matrix["indptr"][:]
        gene_umis = np.zeros(feature_count, dtype=np.int64)
        grid_umis = np.zeros(side * side, dtype=np.int64)
        represented = np.zeros(side * side, dtype=bool)
        total_umis = 0
        for start in range(0, barcode_count, 65536):
            stop = min(start + 65536, barcode_count)
            positions = barcode_grid_indices(matrix["barcodes"][start:stop], side)
            if np.any(represented[positions]) or len(np.unique(positions)) != len(positions):
                raise ValueError(f"Duplicate matrix barcode in {path}")
            first, last = int(pointers[start]), int(pointers[stop])
            values = matrix["data"][first:last].astype(np.int64)
            indices = matrix["indices"][first:last]
            offsets = pointers[start:stop + 1] - first
            cumulative = np.concatenate(([0], np.cumsum(values, dtype=np.int64)))
            per_barcode = np.diff(cumulative[offsets])
            grid_umis[positions] = per_barcode
            represented[positions] = True
            gene_umis += np.bincount(indices, weights=values, minlength=feature_count).astype(np.int64)
            total_umis += int(values.sum(dtype=np.int64))
        if int(grid_umis.sum()) != total_umis or int(gene_umis.sum()) != total_umis:
            raise ValueError(f"Inconsistent gene/barcode sums in {path}")
        return {
            "shape": f"{feature_count} x {barcode_count}",
            "feature_count": feature_count,
            "barcode_count": barcode_count,
            "nonzero_entries": int(matrix["data"].shape[0]),
            "total_umis": total_umis,
            "feature_ids": feature_ids,
            "feature_names": feature_names,
            "software_version": as_text(handle.attrs["software_version"]),
            "grid_umis": grid_umis,
            "represented": represented,
            "gene_umis": gene_umis,
        }


def tissue_positions_stats(path, side):
    mask = np.zeros((side, side), dtype=np.uint8)
    seen = np.zeros(side * side, dtype=bool)
    in_tissue = 0
    for batch in parquet.ParquetFile(path).iter_batches(
        columns=["in_tissue", "array_row", "array_col"], batch_size=250_000
    ):
        values = {name: batch.column(name).to_numpy() for name in batch.schema.names}
        rows = values["array_row"].astype(np.intp)
        columns = values["array_col"].astype(np.intp)
        indices = rows * side + columns
        if np.any(indices >= side * side) or np.any(seen[indices]):
            raise ValueError(f"Invalid or duplicate spatial coordinate in {path}")
        seen[indices] = True
        mask[rows, columns] = values["in_tissue"]
        in_tissue += int(np.sum(values["in_tissue"]))
    if not np.all(seen):
        raise ValueError(f"Incomplete tissue-position grid in {path}")
    return mask, in_tissue


def feature_slice_metadata(path):
    with h5py.File(path, "r") as handle:
        metadata = json.loads(handle.attrs["metadata_json"])
        metadata["software_version"] = as_text(handle.attrs["software_version"])
        return metadata


def correlation(left, right):
    if np.std(left) == 0 or np.std(right) == 0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def compare_matrices(size, kind, local, official):
    common = local["represented"] & official["represented"]
    union = local["represented"] | official["represented"]
    difference = local["grid_umis"] - official["grid_umis"]
    return {
        "bin_size_um": size,
        "matrix_type": kind,
        "features_identical_in_order": local["feature_ids"] == official["feature_ids"] and local["feature_names"] == official["feature_names"],
        "local_barcode_columns": local["barcode_count"],
        "official_barcode_columns": official["barcode_count"],
        "common_barcodes": int(common.sum()),
        "local_only_barcodes": int(np.count_nonzero(local["represented"] & ~official["represented"])),
        "official_only_barcodes": int(np.count_nonzero(official["represented"] & ~local["represented"])),
        "barcode_set_jaccard": float(common.sum() / union.sum()),
        "total_umis_local": int(local["grid_umis"].sum()),
        "total_umis_official": int(official["grid_umis"].sum()),
        "zero_umi_barcodes_local": int(np.count_nonzero(local["represented"] & (local["grid_umis"] == 0))),
        "zero_umi_barcodes_official": int(np.count_nonzero(official["represented"] & (official["grid_umis"] == 0))),
        "local_only_barcodes_with_nonzero_umis": int(np.count_nonzero(local["represented"] & ~official["represented"] & (local["grid_umis"] > 0))),
        "total_umi_difference_percent": float(difference.sum() / official["grid_umis"].sum() * 100),
        "barcode_total_umi_pearson_common": correlation(local["grid_umis"][common], official["grid_umis"][common]),
        "barcode_total_umi_pearson_full_grid": correlation(local["grid_umis"], official["grid_umis"]),
        "different_bin_umi_totals_common": int(np.count_nonzero(difference[common])),
        "sum_absolute_bin_umi_difference_full_grid": int(np.abs(difference).sum()),
        "gene_pseudobulk_umi_pearson": correlation(local["gene_umis"], official["gene_umis"]),
        "gene_pseudobulk_log1p_pearson": correlation(np.log1p(local["gene_umis"]), np.log1p(official["gene_umis"])),
    }


def compare_metrics(metric_sets):
    rows = []
    keys = sorted(set().union(*(set(metrics) for metrics in metric_sets.values())))
    for key in keys:
        local = metric_sets["local_run"].get(key, "")
        official = metric_sets["official_download"].get(key, "")
        if key == "Sample ID":
            rows.append({"metric": key, "unit": "text", "local_run": local, "official_download": official})
            continue
        scale = 100 if key in PERCENT_METRICS else 1
        left = float(local) * scale if local else ""
        right = float(official) * scale if official else ""
        rows.append({
            "metric": key,
            "unit": "percent; difference is percentage points" if scale == 100 else "native",
            "local_run": left,
            "official_download": right,
            "difference_local_minus_official": left - right if local and official else "",
            "relative_difference_percent": (float(local) - float(official)) / float(official) * 100 if local and official and float(official) else "",
        })
    return rows


def write_figures(output_dir, masks):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    figure, axes = plt.subplots(1, 3, figsize=(13, 5), constrained_layout=True)
    rows = []
    for axis, size in zip(axes, SIZES):
        local = masks[("local_run", size)]
        official = masks[("official_download", size)]
        agreement = local + official * 2
        intersection = int(np.count_nonzero(agreement == 3))
        union = int(np.count_nonzero(agreement))
        jaccard = intersection / union if union else 1.0
        rows.append({
            "bin_size_um": size,
            "both_in_tissue": intersection,
            "local_only_in_tissue": int(np.count_nonzero(agreement == 1)),
            "official_only_in_tissue": int(np.count_nonzero(agreement == 2)),
            "in_tissue_jaccard": jaccard,
        })
        axis.imshow(agreement, cmap=ListedColormap(["white", "#e69f00", "#56b4e9", "#444444"]), vmin=0, vmax=3, interpolation="nearest")
        axis.set_title(f"{size} µm tissue mask\nJaccard {jaccard:.4f}")
        axis.set_xlabel("Array column")
        axis.set_ylabel("Array row")
    figure.suptitle("Gray: both; orange: local only; blue: official only")
    figure.savefig(output_dir / "tissue_mask_comparison.png", dpi=180)
    plt.close(figure)
    return rows


def compare_images():
    rows = []
    for name in ("cytassist_image.tiff", "tissue_hires_image.png", "tissue_lowres_image.png"):
        digests = {}
        for label, root in RUNS.items():
            outs = root / "outs" if (root / "outs").is_dir() else root
            path = outs / "spatial" / name
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                    digest.update(chunk)
            digests[label] = digest.hexdigest()
        rows.append({
            "image": name,
            "local_sha256": digests["local_run"],
            "official_sha256": digests["official_download"],
            "byte_identical": digests["local_run"] == digests["official_download"],
        })
    return rows


def write_assessment(path, run_info, metric_sets, matrix_rows, mask_rows, image_rows):
    local_metrics = metric_sets["local_run"]
    official_metrics = metric_sets["official_download"]
    masks = {row["bin_size_um"]: row for row in mask_rows}
    matrices = {(row["bin_size_um"], row["matrix_type"]): row for row in matrix_rows}
    lines = [
        "# Space Ranger 结果比较", "",
        f"比较 `{RUNS['local_run'].relative_to(PROJECT_ROOT)}`（Space Ranger {run_info['local_run']['version']}）与 `{RUNS['official_download'].relative_to(PROJECT_ROOT)}`（Space Ranger {run_info['official_download']['version']}）。",
        f"两份输出报告相同的 `{int(local_metrics['Number of Reads']):,}` 条 reads。",
        "", "## 结论", "",
        "- 两份结果都对应 Mouse Brain Visium HD CytAssist 样本和 slide `H1-8CPZTMY`；mapping/QC 指标非常接近。",
        f"- tissue mask Jaccard：2 µm `{masks[2]['in_tissue_jaccard']:.4f}`，8 µm `{masks[8]['in_tissue_jaccard']:.4f}`，16 µm `{masks[16]['in_tissue_jaccard']:.4f}`。差异可见 `tissue_mask_comparison.png`。",
        f"- 8 µm raw bin UMI 总量相关 `{matrices[(8, 'raw')]['barcode_total_umi_pearson_common']:.6f}`，全基因 pseudobulk UMI 相关 `{matrices[(8, 'raw')]['gene_pseudobulk_umi_pearson']:.9f}`；整体信号高度一致，但逐基因逐 bin 矩阵不是完全相同。",
        f"- 输出中的 CytAssist/组织高低分辨率图像 SHA256 核验：`{sum(row['byte_identical'] for row in image_rows)}/3` 个文件逐字节相同；详见 `image_comparison.csv`。",
        f"- 本地 4.0.1 结果报告 `{local_metrics.get('Number of Cells', 'N/A')}` 个细胞且有 segmentation 文件；官网 3.0.0 输出没有相同口径的细胞级结果。",
        "- 最重要的流程差异：官网下载 HTML 命令包含 `--loupe-alignment`，本地 4.0.1 invocation 未使用手工 alignment；tissue mask 与 bin 数的差异不能直接解释为表达定量失败。",
        "", "## 关键 QC", "",
        "| 指标 | 本地 4.0.1 | 官网 3.0.0 | 本地减官网 |",
        "|---|---:|---:|---:|",
    ]
    for key in ("Valid Barcodes", "Sequencing Saturation", "Reads Mapped Confidently to Probe Set", "Fraction Reads in Squares Under Tissue"):
        left = float(local_metrics[key]) * 100
        right = float(official_metrics[key]) * 100
        lines.append(f"| {key} | {left:.4f}% | {right:.4f}% | {left - right:+.4f} 个百分点 |")
    lines.extend(["", "### Raw bin UMI 与 tissue-bin 数", "", "| bin | raw UMI 4.0.1 | raw UMI 3.0.0 | UMI 差异 | tissue bins 4.0.1 | tissue bins 3.0.0 |", "|---|---:|---:|---:|---:|---:|"])
    for size in SIZES:
        row = matrices[(size, "raw")]
        metric_key = f"Number of {'Squares' if size == 2 else 'Bins'} Under Tissue {size} µm"
        lines.append(f"| {size} µm | {row['total_umis_local']:,} | {row['total_umis_official']:,} | {row['total_umi_difference_percent']:+.4f}% | {int(float(local_metrics[metric_key])):,} | {int(float(official_metrics[metric_key])):,} |")
    lines.extend([
        "", "## 差异原因与限制", "",
        "- 软件版本不同（4.0.1 对 3.0.0）；官网命令使用 Mouse Transcriptome Probe Set v2.0。本地 `outs/probe_set.csv` 文件头同样标为 `Visium Mouse Transcriptome Probe Set v2.0`、mm10、2020-A，因此不是仅因文件名不同就能认定探针版本不同；官网原始探针文件未在项目内，无法证明探针序列逐字节相同。",
        "- 3.0.0 使用 Loupe 手工对齐 JSON，本地 4.0.1 使用自动配准。不同 image registration/tissue detection 是组织 mask 和 tissue-bin 数差异的合理候选原因；未控制变量重跑，不能将各类差异分别定量归因。",
        f"- raw UMI 总量也有 {matrices[(2, 'raw')]['total_umi_difference_percent']:+.4f}% 的差异，不应全部归因于组织检测；barcode/UMI 处理的软件版本变化也是候选原因，但本次未核实版本变更文档。",
        "- 矩阵的 ordered feature ID/名称逐项比较；barcodes 按 array row/column 空间坐标对齐，缺失矩阵列按零计数。bin UMI 和 pseudobulk 相关见 `matrix_comparison.csv`；不声称全量 gene×bin 计数完全一致。",
        f"- 4.0.1 raw 矩阵保留完整空间网格，3.0.0 raw 只包含其中一部分 barcode 列。2 µm 多出的 `{int(matrices[(2, 'raw')]['local_only_barcodes']):,}` 列中，仅 `{int(matrices[(2, 'raw')]['local_only_barcodes_with_nonzero_umis']):,}` 列有非零 UMI，其余是零 UMI 列；不能把 raw 列数增加理解为多出相应数量的组织 bin 或细胞。",
        "- 相同 read 总数和相近 QC 支持二者来自同一公开数据集，但没有原始 FASTQ/reference/probe-set 两边文件的 SHA256，不能声称输入文件字节完全相同。",
        "- 本比较仅读取现有输出，不重跑 Space Ranger、不改动原始结果或 Conda 环境。",
        "", "## 直接复现教程", "",
        "- 用户提供的 10x 页面提供该样本处理后的官方下载结果；`out_10x/Visium_HD_Mouse_Brain_web_summary.html` 内嵌此数据集实际 `spaceranger count` 命令，可直接查看输入参数。",
        "- Space Ranger 官方 `count` 文档解释命令参数。针对这个数据集的运行说明见 `../tutorial.md`；手工运行 3.0.0 命令还需要官方下载的 Loupe alignment JSON，单有输出文件并不能重建该输入。",
        "", "## 结果文件", "",
        "- `metrics_comparison.csv`：QC 指标逐项对比。",
        "- `matrix_comparison.csv`：6 个 resolution/type 组合的矩阵维度、UMI 总数、barcode/gene 对齐统计。",
        "- `tissue_mask_comparison.csv`、`tissue_mask_comparison.png`：tissue 区域交集和图示。",
        "- `run_metadata.csv`：软件版本、原始命令、slide/layout 和告警。", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Compare local and official Space Ranger outputs.")
    parser.add_argument("--results-dir", type=Path, default=SCRIPT_PATH.parents[1] / "results")
    output_dir = parser.parse_args().results_dir.resolve()
    if not output_dir.is_relative_to(SCRIPT_PATH.parents[1]):
        raise ValueError("Results must stay inside spa04/compare")
    output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir))
    metric_sets = {}
    run_info = {}
    matrices = {}
    matrix_rows = []
    masks = {}

    for label, root in RUNS.items():
        outs = root / "outs" if (root / "outs").is_dir() else root
        metric_path = outs / "metrics_summary.csv"
        if not metric_path.exists():
            metric_path = root / "Visium_HD_Mouse_Brain_metrics_summary.csv"
        metric_sets[label] = read_metrics(metric_path)
        summary_path = outs / "web_summary.html"
        if not summary_path.exists():
            summary_path = root / "Visium_HD_Mouse_Brain_web_summary.html"
        summary = read_summary(summary_path)
        feature_path = outs / "feature_slice.h5"
        if not feature_path.exists():
            feature_path = root / "Visium_HD_Mouse_Brain_feature_slice.h5"
        metadata = feature_slice_metadata(feature_path)
        layout = json.loads(metadata["hd_layout_json"])
        run_info[label] = {
            "version": metadata["software_version"],
            "sample": summary["sample"],
            "command": summary["command"],
            "alarm_count": len(summary["alarms"]),
            "slide_uid": layout.get("slide_uid", ""),
            "layout_transform": layout.get("transform", []),
            "registration_transform": metadata["transform_matrices"]["spot_colrow_to_microscope_colrow"],
            "cell_outputs": sorted(path.name for path in (outs / "segmented_outputs").glob("*")) if (outs / "segmented_outputs").is_dir() else [],
            "summary_path": str(summary_path.relative_to(PROJECT_ROOT)),
        }
        for size in SIZES:
            side = (6700 + size - 1) // size
            matrix_folder = outs / "binned_outputs" / f"square_{size:03d}um"
            spatial_path = matrix_folder / "spatial" / "tissue_positions.parquet"
            tissue_mask, positions_in_tissue = tissue_positions_stats(spatial_path, side)
            masks[(label, size)] = tissue_mask
            for kind in ("raw", "filtered"):
                matrix_path = matrix_folder / f"{kind}_feature_bc_matrix.h5"
                stats = matrix_stats(matrix_path, side)
                matrices[(label, size, kind)] = stats
                matrix_rows.append({
                    "run": label,
                    "space_ranger_version": stats["software_version"],
                    "bin_size_um": size,
                    "matrix_type": kind,
                    "matrix_shape_features_x_barcodes": stats["shape"],
                    "feature_count": stats["feature_count"],
                    "barcode_count": stats["barcode_count"],
                    "nonzero_entries": stats["nonzero_entries"],
                    "total_umis": stats["total_umis"],
                    "positions_in_tissue": positions_in_tissue,
                })

    comparison_rows = []
    for size in SIZES:
        for kind in ("raw", "filtered"):
            local = matrices[("local_run", size, kind)]
            official = matrices[("official_download", size, kind)]
            if local["feature_ids"] != official["feature_ids"] or local["feature_names"] != official["feature_names"]:
                raise ValueError(f"Feature definitions differ at {size} µm {kind} matrix")
            comparison_rows.append(compare_matrices(size, kind, local, official))

    write_csv(output_dir / "metrics_comparison.csv", compare_metrics(metric_sets))
    write_csv(output_dir / "matrix_comparison.csv", comparison_rows)
    mask_rows = write_figures(output_dir, masks)
    write_csv(output_dir / "tissue_mask_comparison.csv", mask_rows)
    write_csv(output_dir / "matrix_inventory.csv", matrix_rows)
    image_rows = compare_images()
    write_csv(output_dir / "image_comparison.csv", image_rows)

    metadata_rows = [{
        "run": label,
        "space_ranger_version": info["version"],
        "sample_id": info["sample"].get("id", ""),
        "dataset_description": info["sample"].get("description", ""),
        "command_from_web_summary": info["command"],
        "web_summary_alarm_count": info["alarm_count"],
        "cell_segmentation_outputs": "; ".join(info["cell_outputs"]),
        "slide_uid": info["slide_uid"],
        "automatic_layout_transform": json.dumps(info["layout_transform"]),
        "spot_to_microscope_transform": json.dumps(info["registration_transform"]),
        "web_summary_path": info["summary_path"],
    } for label, info in run_info.items()]
    write_csv(output_dir / "run_metadata.csv", metadata_rows)
    write_assessment(output_dir / "assessment.md", run_info, metric_sets, comparison_rows, mask_rows, image_rows)
    print(f"Wrote Space Ranger comparison to {output_dir}")


if __name__ == "__main__":
    main()
