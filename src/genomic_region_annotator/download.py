from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import requests


@dataclass(frozen=True)
class EnsemblGTF:
    release: int
    species: str = "homo_sapiens"
    assembly: str | None = None  # optional filter for filename selection

    @property
    def base_dir_url(self) -> str:
        return f"https://ftp.ensembl.org/pub/release-{self.release}/gtf/{self.species}/"


def _run(cmd: list[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stdout)


def _discover_gtf_filename(release: int, species: str, assembly: str | None = None) -> str:
    """
    Discover the appropriate *.gtf.gz filename for a species and Ensembl release
    by listing the Ensembl HTTPS directory and selecting a best candidate.

    If assembly is provided, it is used as a case-insensitive substring filter.
    """
    url = f"https://ftp.ensembl.org/pub/release-{release}/gtf/{species}/"
    r = requests.get(url, timeout=60)
    r.raise_for_status()

    # Directory listing contains href="...gtf.gz"
    files = sorted(set(re.findall(r'href="([^"]+\.gtf\.gz)"', r.text)))
    if not files:
        raise RuntimeError(f"No .gtf.gz files found at {url}")

    if assembly:
        cand = [f for f in files if assembly.lower() in f.lower()]
        if not cand:
            raise RuntimeError(f"No .gtf.gz matched assembly='{assembly}' at {url}")
        # Heuristic: shortest is usually the canonical one
        return sorted(cand, key=len)[0]

    # No assembly requested: choose a canonical default.
    # Heuristic: shortest filename is often the primary annotation file.
    return sorted(files, key=len)[0]


def download_gtf(
    release: int,
    out_dir: str,
    species: str = "homo_sapiens",
    assembly: str | None = None,
) -> Path:
    """
    Download Ensembl GTF for the given release/species via HTTPS.
    Returns path to uncompressed .gtf file.

    - species: Ensembl species directory name (e.g. homo_sapiens, mus_musculus)
    - assembly: optional substring filter to select among multiple GTFs
    """
    spec = EnsemblGTF(release=release, species=species, assembly=assembly)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    filename = _discover_gtf_filename(release=release, species=species, assembly=assembly)
    https_url = spec.base_dir_url + filename

    gz_path = out / filename
    gtf_path = out / filename.replace(".gz", "")

    # Reuse if already downloaded/unzipped
    if gtf_path.exists():
        return gtf_path

    # Download (prefer curl, fallback to wget)
    if shutil.which("curl"):
        _run(["curl", "-L", "--retry", "5", "--retry-delay", "3", "-o", str(gz_path), https_url])
    elif shutil.which("wget"):
        _run(["wget", "-c", "-O", str(gz_path), https_url])
    else:
        raise RuntimeError("Neither 'curl' nor 'wget' found on this system.")

    # Unzip
    _run(["gunzip", "-f", str(gz_path)])

    if not gtf_path.exists():
        raise RuntimeError(f"Expected {gtf_path} after gunzip, but it was not found.")

    return gtf_path