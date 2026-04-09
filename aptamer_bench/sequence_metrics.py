"""
Per-sequence quality metrics for aptamer evaluation.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from aptamer_bench.sequence import AptamerSequence


@dataclass
class SequenceMetricsResult:
    """Container for all computed per-sequence metrics."""
    aptamer_id: str
    sequence: str
    length: int
    gc_content: float
    linguistic_complexity: float
    has_g_quadruplex: bool
    nucleotide_frequencies: dict[str, float]
    kd_nm: float | None = None

    def to_dict(self) -> dict:
        return {
            "aptamer_id": self.aptamer_id,
            "sequence": self.sequence,
            "length": self.length,
            "gc_content": round(self.gc_content, 4),
            "linguistic_complexity": round(self.linguistic_complexity, 4),
            "has_g_quadruplex": self.has_g_quadruplex,
            "kd_nm": self.kd_nm,
            **{f"freq_{k}": round(v, 4) for k, v in self.nucleotide_frequencies.items()},
        }


class SequenceMetrics:
    """
    Compute per-sequence quality metrics for one or more aptamers.

    Usage
    -----
    >>> from aptamer_bench import AptamerSequence, SequenceMetrics
    >>> apt = AptamerSequence("GCCTGTTGTGAGCCTCCTGTCGAA", seq_type="DNA", target="VEGF")
    >>> sm = SequenceMetrics()
    >>> result = sm.compute(apt)
    >>> print(result.gc_content)
    """

    def compute(self, aptamer: AptamerSequence) -> SequenceMetricsResult:
        """
        Compute all metrics for a single aptamer.

        Parameters
        ----------
        aptamer : AptamerSequence

        Returns
        -------
        SequenceMetricsResult
        """
        return SequenceMetricsResult(
            aptamer_id=aptamer.aptamer_id,
            sequence=aptamer.sequence,
            length=aptamer.length,
            gc_content=aptamer.gc_content,
            linguistic_complexity=self._linguistic_complexity(aptamer.sequence),
            has_g_quadruplex=aptamer.has_g_quadruplex_motif(),
            nucleotide_frequencies=aptamer.nucleotide_frequencies,
            kd_nm=aptamer.kd_nm,
        )

    def compute_batch(self, aptamers: list[AptamerSequence]) -> list[SequenceMetricsResult]:
        """Compute metrics for a list of aptamers."""
        return [self.compute(a) for a in aptamers]

    @staticmethod
    def _linguistic_complexity(sequence: str) -> float:
        """
        Compute linguistic complexity (Ziv-Merhav) of a sequence.

        Linguistic complexity measures sequence diversity as the ratio of
        observed unique substrings to the theoretical maximum. Values close
        to 1.0 indicate high complexity (good for aptamers — avoids repetition).

        Parameters
        ----------
        sequence : str
            Nucleotide sequence.

        Returns
        -------
        float
            Complexity score in [0, 1].
        """
        n = len(sequence)
        if n == 0:
            return 0.0

        observed = set()
        for length in range(1, n + 1):
            for i in range(n - length + 1):
                observed.add(sequence[i : i + length])

        # Theoretical maximum number of unique substrings for alphabet size 4
        alphabet_size = 4
        max_substrings = sum(
            min(alphabet_size ** k, n - k + 1) for k in range(1, n + 1)
        )

        return len(observed) / max_substrings if max_substrings > 0 else 0.0

    @staticmethod
    def shannon_entropy(sequence: str) -> float:
        """
        Compute Shannon entropy of a nucleotide sequence.

        Higher entropy indicates more uniform nucleotide distribution.

        Parameters
        ----------
        sequence : str

        Returns
        -------
        float
            Entropy in bits.
        """
        n = len(sequence)
        if n == 0:
            return 0.0
        freqs = [sequence.count(b) / n for b in set(sequence) if sequence.count(b) > 0]
        return -sum(f * math.log2(f) for f in freqs)
