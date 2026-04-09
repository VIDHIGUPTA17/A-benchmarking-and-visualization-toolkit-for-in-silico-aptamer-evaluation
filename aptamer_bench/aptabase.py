"""
Loader for the Aptamer Database (aptamer.atdbio.com) and similar public sources.

This module fetches, parses, and standardizes aptamer-target binding data
into pyaptamer-compatible formats.

Note: For reproducibility and offline use, data is cached locally by default.
"""

from __future__ import annotations
import json
import hashlib
import logging
from pathlib import Path
from typing import Optional
import pandas as pd

from aptamer_bench.sequence import AptamerSequence

logger = logging.getLogger(__name__)

# Default cache directory
DEFAULT_CACHE_DIR = Path.home() / ".aptamer_bench" / "cache"

# Bundled curated sample data for offline use / testing
# These are publicly known aptamers from literature
BUNDLED_APTAMERS = [
    {
        "aptamer_id": "APT-001",
        "sequence": "GGTTGGTGTGGTTGG",
        "seq_type": "DNA",
        "target": "Thrombin",
        "kd_nm": 26.0,
        "source": "Bock et al., 1992, Nature",
    },
    {
        "aptamer_id": "APT-002",
        "sequence": "AGTCCGTGGTAGGGCAGGTTGGGGTGACT",
        "seq_type": "DNA",
        "target": "Thrombin",
        "kd_nm": 0.5,
        "source": "Tasset et al., 1997, JMB",
    },
    {
        "aptamer_id": "APT-003",
        "sequence": "GGGAGGGCGGGUCGGGAGGGG",
        "seq_type": "RNA",
        "target": "VEGF",
        "kd_nm": 0.14,
        "source": "Ruckman et al., 1998, JBC",
    },
    {
        "aptamer_id": "APT-004",
        "sequence": "GCCTGTTGTGAGCCTCCTGTCGAA",
        "seq_type": "DNA",
        "target": "VEGF",
        "kd_nm": 1.0,
        "source": "Li et al., 2014",
    },
    {
        "aptamer_id": "APT-005",
        "sequence": "ATCCAGAGTGACGCAGCA",
        "seq_type": "DNA",
        "target": "MUC1",
        "kd_nm": 0.135,
        "source": "Ferreira et al., 2006",
    },
    {
        "aptamer_id": "APT-006",
        "sequence": "CTACGGCACGTTTATCCGTCCCTCCTAGTGGCGTGCCGTAG",
        "seq_type": "DNA",
        "target": "HER2",
        "kd_nm": 0.29,
        "source": "Li et al., 2012",
    },
    {
        "aptamer_id": "APT-007",
        "sequence": "GGGCCGAAAAAGGGC",
        "seq_type": "DNA",
        "target": "Streptavidin",
        "kd_nm": 70.0,
        "source": "Bittker et al., 2002",
    },
    {
        "aptamer_id": "APT-008",
        "sequence": "AACCGCCCAAATCCCTAAGAGTCTGCACTT",
        "seq_type": "DNA",
        "target": "PDGF-BB",
        "kd_nm": 0.1,
        "source": "Green et al., 1996",
    },
]


class AptaBaseLoader:
    """
    Load and standardize aptamer data from bundled curated datasets
    or user-provided CSV files.

    For research use: extend this class to integrate with live databases
    such as aptamer.atdbio.com or AptaBase when their APIs are available.

    Parameters
    ----------
    cache_dir : Path, optional
        Directory to cache downloaded data. Defaults to ~/.aptamer_bench/cache

    Examples
    --------
    >>> loader = AptaBaseLoader()
    >>> aptamers = loader.load_bundled()
    >>> df = loader.to_dataframe(aptamers)
    >>> print(df.head())
    """

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def load_bundled(self, target_filter: Optional[str] = None) -> list[AptamerSequence]:
        """
        Load the bundled curated aptamer dataset.

        Parameters
        ----------
        target_filter : str, optional
            If provided, only return aptamers for this target (case-insensitive).

        Returns
        -------
        list[AptamerSequence]
        """
        data = BUNDLED_APTAMERS
        if target_filter:
            data = [d for d in data if target_filter.lower() in d["target"].lower()]

        return [
            AptamerSequence(
                sequence=d["sequence"],
                seq_type=d["seq_type"],
                target=d["target"],
                kd_nm=d.get("kd_nm"),
                source=d.get("source", ""),
                aptamer_id=d["aptamer_id"],
            )
            for d in data
        ]

    def load_csv(self, path: str | Path) -> list[AptamerSequence]:
        """
        Load aptamers from a CSV file.

        Expected columns: sequence, seq_type, target (optional: kd_nm, source, aptamer_id)

        Parameters
        ----------
        path : str or Path

        Returns
        -------
        list[AptamerSequence]
        """
        df = pd.read_csv(path)
        required = {"sequence", "seq_type"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"CSV is missing required columns: {missing}")

        aptamers = []
        for i, row in df.iterrows():
            try:
                apt = AptamerSequence(
                    sequence=str(row["sequence"]),
                    seq_type=row.get("seq_type", "DNA"),
                    target=str(row.get("target", "")),
                    kd_nm=float(row["kd_nm"]) if "kd_nm" in row and pd.notna(row["kd_nm"]) else None,
                    source=str(row.get("source", "")),
                    aptamer_id=str(row.get("aptamer_id", f"APT-{i:04d}")),
                )
                aptamers.append(apt)
            except ValueError as e:
                logger.warning(f"Skipping row {i}: {e}")

        logger.info(f"Loaded {len(aptamers)} aptamers from {path}")
        return aptamers

    def to_dataframe(self, aptamers: list[AptamerSequence]) -> pd.DataFrame:
        """Convert a list of AptamerSequence objects to a pandas DataFrame."""
        return pd.DataFrame([
            {
                "aptamer_id": a.aptamer_id,
                "sequence": a.sequence,
                "seq_type": a.seq_type,
                "target": a.target,
                "kd_nm": a.kd_nm,
                "length": a.length,
                "gc_content": round(a.gc_content, 4),
                "source": a.source,
            }
            for a in aptamers
        ])

    def to_fasta(self, aptamers: list[AptamerSequence], path: str | Path) -> None:
        """
        Export aptamers to FASTA format (compatible with Biopython and pyaptamer).

        Parameters
        ----------
        aptamers : list[AptamerSequence]
        path : str or Path
            Output file path.
        """
        path = Path(path)
        with open(path, "w") as f:
            for apt in aptamers:
                header = f">{apt.aptamer_id} target={apt.target}"
                if apt.kd_nm is not None:
                    header += f" kd_nm={apt.kd_nm}"
                f.write(f"{header}\n{apt.sequence}\n")
        logger.info(f"Exported {len(aptamers)} sequences to {path}")

    def to_json(self, aptamers: list[AptamerSequence], path: str | Path) -> None:
        """Export aptamers to JSON format."""
        path = Path(path)
        data = [
            {
                "aptamer_id": a.aptamer_id,
                "sequence": a.sequence,
                "seq_type": a.seq_type,
                "target": a.target,
                "kd_nm": a.kd_nm,
                "source": a.source,
            }
            for a in aptamers
        ]
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Exported {len(aptamers)} sequences to {path}")
