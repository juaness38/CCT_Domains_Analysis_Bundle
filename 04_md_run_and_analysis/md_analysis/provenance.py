"""CCT-OSR1 Analysis Provenance Utilities.

Goal: make each analysis step auditable by emitting a small JSON artifact capturing:
- inputs (paths + size/mtime; optional sha256)
- parameters
- outputs
- code identity (file hashes + best-effort git commit)
- runtime environment (python/platform + key package versions)

Design principles:
- Stdlib-only (no extra deps).
- Best-effort: never fail the analysis because provenance cannot be collected.
- Cheap by default: hashes are computed for small files; large files record metadata unless forced.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass
class FileRecord:
    path: str
    exists: bool
    size_bytes: Optional[int] = None
    mtime_iso: Optional[str] = None
    sha256: Optional[str] = None


@dataclass
class ProvenanceRecord:
    schema_version: str
    kind: str  # "run" | "analysis"
    created_utc: str
    status: str  # "ok" | "error"
    name: str

    argv: List[str] = field(default_factory=list)
    cwd: Optional[str] = None

    inputs: List[FileRecord] = field(default_factory=list)
    outputs: List[FileRecord] = field(default_factory=list)

    params: Dict[str, Any] = field(default_factory=dict)

    code_files: List[FileRecord] = field(default_factory=list)
    git_commit: Optional[str] = None

    python: Dict[str, Any] = field(default_factory=dict)
    packages: Dict[str, Optional[str]] = field(default_factory=dict)
    env: Dict[str, str] = field(default_factory=dict)

    error: Optional[Dict[str, Any]] = None


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _safe_stat(path: Path) -> Tuple[bool, Optional[int], Optional[str]]:
    try:
        st = path.stat()
        mtime = datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(timespec="seconds")
        return True, int(st.st_size), mtime
    except FileNotFoundError:
        return False, None, None
    except Exception:
        return False, None, None


def _sha256_file(path: Path, *, max_bytes: Optional[int] = None) -> Optional[str]:
    try:
        h = hashlib.sha256()
        with path.open("rb") as f:
            if max_bytes is None:
                for chunk in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(chunk)
            else:
                remaining = max_bytes
                while remaining > 0:
                    chunk = f.read(min(1024 * 1024, remaining))
                    if not chunk:
                        break
                    h.update(chunk)
                    remaining -= len(chunk)
        return h.hexdigest()
    except Exception:
        return None


def file_record(
    path: os.PathLike[str] | str,
    *,
    hash_if_small_bytes: int = 50 * 1024 * 1024,
    force_hash: bool = False,
    max_hash_bytes: Optional[int] = None,
) -> FileRecord:
    p = Path(path)
    exists, size, mtime = _safe_stat(p)

    sha = None
    if exists and size is not None:
        if force_hash or size <= hash_if_small_bytes:
            sha = _sha256_file(p, max_bytes=max_hash_bytes)

    return FileRecord(
        path=str(p.resolve()) if exists else str(p),
        exists=exists,
        size_bytes=size,
        mtime_iso=mtime,
        sha256=sha,
    )


def best_effort_git_commit(repo_root: os.PathLike[str] | str) -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            commit = result.stdout.strip()
            return commit or None
    except Exception:
        pass
    return None


def get_package_versions(names: Sequence[str]) -> Dict[str, Optional[str]]:
    versions: Dict[str, Optional[str]] = {n: None for n in names}
    try:
        from importlib.metadata import version as pkg_version  # py>=3.8

        for n in names:
            try:
                versions[n] = pkg_version(n)
            except Exception:
                versions[n] = None
    except Exception:
        # Very old python; keep None.
        pass
    return versions


def write_provenance_json(
    out_path: os.PathLike[str] | str,
    record: ProvenanceRecord,
) -> None:
    try:
        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(asdict(record), indent=2, sort_keys=True), encoding="utf-8")
    except Exception:
        # Never raise from provenance.
        return


def make_run_provenance(
    *,
    name: str,
    repo_root: os.PathLike[str] | str,
    inputs: Sequence[os.PathLike[str] | str] = (),
    params: Optional[Dict[str, Any]] = None,
    env_prefixes: Sequence[str] = ("CCT_", "MDIFP_"),
    extra_packages: Sequence[str] = (),
    status: str = "ok",
    error: Optional[Dict[str, Any]] = None,
) -> ProvenanceRecord:
    packages = [
        "numpy",
        "pandas",
        "mdtraj",
        "MDAnalysis",
        "scipy",
        "plotly",
        "kaleido",
        "scikit-learn",
    ] + list(extra_packages)

    env: Dict[str, str] = {}
    for k, v in os.environ.items():
        if any(k.startswith(pfx) for pfx in env_prefixes):
            env[k] = v

    rec = ProvenanceRecord(
        schema_version="1.0",
        kind="run",
        created_utc=_utc_now_iso(),
        status=status,
        name=name,
        argv=list(sys.argv),
        cwd=str(Path.cwd()),
        inputs=[file_record(p) for p in inputs],
        outputs=[],
        params=params or {},
        code_files=[],
        git_commit=best_effort_git_commit(repo_root),
        python={
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
        packages=get_package_versions(packages),
        env=env,
        error=error,
    )
    return rec


def make_analysis_provenance(
    *,
    name: str,
    repo_root: os.PathLike[str] | str,
    inputs: Sequence[os.PathLike[str] | str] = (),
    outputs: Sequence[os.PathLike[str] | str] = (),
    params: Optional[Dict[str, Any]] = None,
    code_paths: Sequence[os.PathLike[str] | str] = (),
    status: str = "ok",
    error: Optional[Dict[str, Any]] = None,
) -> ProvenanceRecord:
    rec = ProvenanceRecord(
        schema_version="1.0",
        kind="analysis",
        created_utc=_utc_now_iso(),
        status=status,
        name=name,
        argv=list(sys.argv),
        cwd=str(Path.cwd()),
        inputs=[file_record(p) for p in inputs],
        outputs=[file_record(p, hash_if_small_bytes=50 * 1024 * 1024) for p in outputs],
        params=params or {},
        code_files=[file_record(p, force_hash=True) for p in code_paths],
        git_commit=best_effort_git_commit(repo_root),
        python={
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
        packages=get_package_versions([
            "numpy",
            "pandas",
            "mdtraj",
            "MDAnalysis",
            "scipy",
            "plotly",
        ]),
        env={},
        error=error,
    )
    return rec
