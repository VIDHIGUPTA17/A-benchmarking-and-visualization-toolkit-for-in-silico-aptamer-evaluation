"""
Standardized benchmark datasets for comparing in-silico aptamer generation algorithms.

BenchmarkDataset wraps a set of reference aptamers and provides consistent
evaluation against algorithm-generated candidates — directly addressing the
gap identified in the pyaptamer roadmap.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import pandas as pd

from aptamer_bench.sequence import AptamerSequence
from aptamer_bench.sequence_metrics import SequenceMetrics
from diversity_metrics import DiversityMetrics


@dataclass
class BenchmarkResult:
    """Results of benchmarking a set of generated aptamers against a reference."""
    algorithm_name: str
    n_generated: int
    diversity_summary: dict
    per_sequence_metrics: list[dict]
    reference_kd_mean_nm: Optional[float] = None
    notes: str = ""

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.per_sequence_metrics)

    def __repr__(self) -> str:
        return (
            f"BenchmarkResult(algorithm='{self.algorithm_name}', "
            f"n={self.n_generated}, "
            f"redundancy={self.diversity_summary.get('redundancy_rate', 'N/A'):.3f}, "
            f"mean_distance={self.diversity_summary.get('mean_pairwise_distance_normalized', 'N/A'):.3f})"
        )


class BenchmarkDataset:
    """
    Wraps reference aptamer data and evaluates algorithm-generated candidates.

    This is the core benchmarking class designed to be compatible with
    pyaptamer algorithm outputs.

    Parameters
    ----------
    reference_aptamers : list[AptamerSequence]
        Known, validated aptamers to serve as ground truth / comparison baseline.
    name : str
        Human-readable name for this benchmark (e.g., "Thrombin-Binding-DNA").

    Examples
    --------
    >>> from aptamer_bench.datasets import AptaBaseLoader, BenchmarkDataset
    >>> loader = AptaBaseLoader()
    >>> refs = loader.load_bundled(target_filter="Thrombin")
    >>> bench = BenchmarkDataset(reference_aptamers=refs, name="Thrombin")
    >>> # Suppose your algorithm generated these candidates:
    >>> candidates = [AptamerSequence("GGTTGGTGTGGTTGG", "DNA")]
    >>> result = bench.evaluate(candidates, algorithm_name="MyAlgo")
    >>> print(result)
    """

    def __init__(
        self,
        reference_aptamers: list[AptamerSequence],
        name: str = "unnamed",
    ):
        self.reference_aptamers = reference_aptamers
        self.name = name
        self._seq_metrics = SequenceMetrics()
        self._div_metrics = DiversityMetrics()

    @property
    def n_reference(self) -> int:
        return len(self.reference_aptamers)

    def reference_summary(self) -> dict:
        """Summary statistics for the reference dataset."""
        return {
            "name": self.name,
            "n_reference": self.n_reference,
            "diversity": self._div_metrics.summary(self.reference_aptamers),
        }

    def evaluate(
        self,
        generated: list[AptamerSequence],
        algorithm_name: str = "unknown",
        notes: str = "",
    ) -> BenchmarkResult:
        """
        Evaluate a set of algorithm-generated aptamers.

        Computes:
        - Per-sequence metrics (GC content, complexity, G-quadruplex)
        - Population diversity (redundancy, pairwise distance, length/GC distributions)

        Parameters
        ----------
        generated : list[AptamerSequence]
            Candidate aptamers produced by an in-silico algorithm.
        algorithm_name : str
            Label for the algorithm being evaluated.
        notes : str
            Free-text notes about the run.

        Returns
        -------
        BenchmarkResult
        """
        seq_results = self._seq_metrics.compute_batch(generated)
        div_summary = self._div_metrics.summary(generated)

        ref_kds = [a.kd_nm for a in self.reference_aptamers if a.kd_nm is not None]
        ref_kd_mean = sum(ref_kds) / len(ref_kds) if ref_kds else None

        return BenchmarkResult(
            algorithm_name=algorithm_name,
            n_generated=len(generated),
            diversity_summary=div_summary,
            per_sequence_metrics=[r.to_dict() for r in seq_results],
            reference_kd_mean_nm=ref_kd_mean,
            notes=notes,
        )

    def compare_algorithms(
        self,
        algorithm_outputs: dict[str, list[AptamerSequence]],
    ) -> pd.DataFrame:
        """
        Compare multiple algorithms side by side.

        Parameters
        ----------
        algorithm_outputs : dict[str, list[AptamerSequence]]
            Keys are algorithm names, values are lists of generated aptamers.

        Returns
        -------
        pd.DataFrame
            One row per algorithm with summary metrics.
        """
        rows = []
        for algo_name, sequences in algorithm_outputs.items():
            result = self.evaluate(sequences, algorithm_name=algo_name)
            div = result.diversity_summary
            rows.append({
                "algorithm": algo_name,
                "n_generated": result.n_generated,
                "n_unique": div["n_unique"],
                "redundancy_rate": round(div["redundancy_rate"], 4),
                "mean_pairwise_distance": round(
                    div["mean_pairwise_distance_normalized"], 4
                ),
                "gc_content_mean": round(div["gc_content"]["mean"], 4),
                "gc_content_std": round(div["gc_content"]["std"], 4),
                "length_mean": round(div["length"]["mean"], 2),
            })
        return pd.DataFrame(rows).sort_values("redundancy_rate")
