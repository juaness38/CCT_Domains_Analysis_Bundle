from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Sequence

import pandas as pd

import _replot_md_controls_rmsd_rmsf as base


OUTPUT_ROOT = base.OUTPUT_ROOT / "rmsf_only_v2"


def _dashboard_html(records: Sequence[Dict[str, Any]], output_html: Path) -> None:
    cards = []
    for record in records:
        rel_html = Path(record["rmsf_html_name"]).name
        rel_png = Path(record["rmsf_png_name"]).name
        cards.append(
            f"""
            <article>
              <h2>{record['system_title']}</h2>
              <p><strong>Input RMSF CSV:</strong> {record['rmsf_csv']}</p>
              <p><strong>RMSF mode:</strong> {record['rmsf_mode']}</p>
              <p><strong>RMSF source:</strong> {record['rmsf_runtime_source']}</p>
              <p><strong>Sampling stride:</strong> {record['rmsf_sampling_stride']}</p>
              <p><strong>Sampled frames:</strong> {record['rmsf_sampled_frames']}</p>
              <p><a href=\"{rel_html}\">Interactive RMSF HTML</a></p>
              <figure>
                <img src=\"{rel_png}\" alt=\"{record['system_title']} RMSF\" />
                <figcaption>Corrected self-aligned C-alpha RMSF</figcaption>
              </figure>
            </article>
            """
        )

    html = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>Control MD RMSF Only Dashboard</title>
  <style>
    body {{ font-family: Arial, sans-serif; background: linear-gradient(180deg, #eef6fb 0%, #f8fbfd 100%); color: #14364c; margin: 0; }}
    main {{ max-width: 1800px; margin: 0 auto; padding: 28px 24px 40px; }}
    h1 {{ margin: 0 0 10px; font-size: 38px; }}
    p {{ line-height: 1.55; }}
    article {{ background: white; border-radius: 20px; padding: 22px; box-shadow: 0 12px 34px rgba(20, 54, 76, 0.08); margin-top: 22px; }}
    figure {{ margin: 16px 0 0; background: #f7fafc; border-radius: 14px; padding: 12px; }}
    img {{ width: 100%; height: auto; border-radius: 10px; display: block; }}
    figcaption {{ padding-top: 10px; font-size: 14px; }}
    a {{ color: #1e6ea1; }}
  </style>
</head>
<body>
  <main>
    <h1>Unified MD RMSF Only Dashboard</h1>
    <p>Versioned RMSF-only export. All figures are regenerated from raw trajectories using self-aligned C-alpha RMSF and are isolated from the mixed RMSD/RMSF bundle.</p>
    {''.join(cards)}
  </main>
</body>
</html>
"""
    output_html.write_text(html, encoding="utf-8")


def _process_replica(spec: base.ReplicaSpec, output_dir: Path) -> Dict[str, Any]:
    rmsf_csv = base._find_csv(spec.source_root, "rmsf", spec.preferred_prefix)
    rmsf_df, rmsf_mode, rmsf_meta = base._build_corrected_rmsf_dataframe(spec)

    system_title = base._system_title(spec.system_label, spec.replica_label)
    base_name = base._slugify(system_title)

    rmsf_fig = base._plot_rmsf(system_title, rmsf_df, rmsf_mode)

    rmsf_png = output_dir / f"{base_name}_rmsf_ca.png"
    rmsf_html = output_dir / f"{base_name}_rmsf_ca.html"
    rmsf_data_csv = output_dir / f"{base_name}_rmsf_ca_corrected.csv"

    rmsf_fig.write_html(str(rmsf_html))
    base._write_rmsf_png(system_title, rmsf_df, rmsf_png)
    rmsf_df.to_csv(rmsf_data_csv, index=False)

    return {
        "system_label": spec.system_label,
        "replica_label": spec.replica_label,
        "system_title": system_title,
        "source_root": str(spec.source_root),
        "rmsf_csv": str(rmsf_csv),
        "rmsf_mode": rmsf_mode,
        "rmsf_runtime_source": "raw_dcd_self_aligned_ca",
        "rmsf_sampling_stride": rmsf_meta["rmsf_sampling_stride"],
        "rmsf_sampled_frames": rmsf_meta["rmsf_sampled_frames"],
        "runtime_source": spec.runtime_source,
        "rmsf_html_name": rmsf_html.name,
        "rmsf_png_name": rmsf_png.name,
        "rmsf_data_csv_name": rmsf_data_csv.name,
    }


def _write_summary(records: Sequence[Dict[str, Any]], output_dir: Path) -> None:
    summary_df = pd.DataFrame(records)
    summary_df.to_csv(output_dir / "replot_manifest.csv", index=False)

    (output_dir / "replot_manifest.json").write_text(
        json.dumps(records, indent=2),
        encoding="utf-8",
    )

    lines = [
        "# Unified MD RMSF Only Replot Manifest",
        "",
        "This document records each resolved RMSF source and the regenerated RMSF-only artifacts.",
        "",
    ]
    for record in records:
        lines.extend(
            [
                f"## {record['system_title']}",
                "",
                f"- Source root: {record['source_root']}",
                f"- Input RMSF CSV: {record['rmsf_csv']}",
                f"- RMSF mode: {record['rmsf_mode']}",
                f"- RMSF source: {record['rmsf_runtime_source']}",
                f"- RMSF sampling stride: {record['rmsf_sampling_stride']}",
                f"- RMSF sampled frames: {record['rmsf_sampled_frames']}",
                f"- Runtime source reference: {record['runtime_source']}",
                f"- RMSF HTML: {record['rmsf_html_name']}",
                f"- RMSF PNG: {record['rmsf_png_name']}",
                f"- Corrected RMSF CSV: {record['rmsf_data_csv_name']}",
                "",
            ]
        )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    records = [_process_replica(spec, OUTPUT_ROOT) for spec in base.REPLICA_SPECS]
    _write_summary(records, OUTPUT_ROOT)
    _dashboard_html(records, OUTPUT_ROOT / "all_systems_rmsf_dashboard.html")
    print(f"[OK] Wrote RMSF-only replots to: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()