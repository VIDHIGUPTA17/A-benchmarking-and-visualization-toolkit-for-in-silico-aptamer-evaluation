"""
Core sequence representation for aptamers.
Supports both DNA and RNA aptamer sequences with validation and basic analysis.
"""

from __future__ import annotations
import re
from dataclasses import dataclass, field
from typing import Literal

# Valid nucleotide alphabets
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")

SeqType = Literal["DNA", "RNA"]


@dataclass
class AptamerSequence:
    """
    Represents a single aptamer sequence with metadata.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence (uppercase). DNA: ACGT, RNA: ACGU.
    seq_type : str
        Either 'DNA' or 'RNA'.
    target : str, optional
        Name/identifier of the binding target (e.g., protein name).
    kd_nm : float, optional
        Dissociation constant (Kd) in nanomolar. Lower = stronger binding.
    source : str, optional
        Database or paper source of this sequence.
    aptamer_id : str, optional
        Unique identifier for this aptamer.

    Examples
    --------
    >>> apt = AptamerSequence("GCCTGTTGTGAGCCTCCTGTCGAA", seq_type="DNA", target="VEGF")
    >>> apt.gc_content
    0.5416666666666666
    >>> apt.length
    24
    """

    sequence: str
    seq_type: SeqType = "DNA"
    target: str = ""
    kd_nm: float | None = None
    source: str = ""
    aptamer_id: str = ""
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.sequence = self.sequence.upper().strip()
        self._validate()

    def _validate(self):
        """Validate sequence characters against the declared alphabet."""
        alphabet = DNA_ALPHABET if self.seq_type == "DNA" else RNA_ALPHABET
        invalid = set(self.sequence) - alphabet
        if invalid:
            raise ValueError(
                f"Invalid characters {invalid} for {self.seq_type} sequence. "
                f"Allowed: {alphabet}"
            )

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------

    @property
    def length(self) -> int:
        """Return length of the sequence."""
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        """Return GC content as a fraction (0.0 to 1.0)."""
        if self.length == 0:
            return 0.0
        gc = self.sequence.count("G") + self.sequence.count("C")
        return gc / self.length

    @property
    def nucleotide_frequencies(self) -> dict[str, float]:
        """Return per-nucleotide frequencies as fractions."""
        bases = "ACGT" if self.seq_type == "DNA" else "ACGU"
        return {b: self.sequence.count(b) / self.length for b in bases}

    # ------------------------------------------------------------------
    # Conversions
    # ------------------------------------------------------------------

    def to_rna(self) -> "AptamerSequence":
        """Convert a DNA aptamer to its RNA equivalent (T → U)."""
        if self.seq_type == "RNA":
            return self
        rna_seq = self.sequence.replace("T", "U")
        return AptamerSequence(
            sequence=rna_seq,
            seq_type="RNA",
            target=self.target,
            kd_nm=self.kd_nm,
            source=self.source,
            aptamer_id=self.aptamer_id,
            metadata=self.metadata,
        )

    def to_dna(self) -> "AptamerSequence":
        """Convert an RNA aptamer to its DNA equivalent (U → T)."""
        if self.seq_type == "DNA":
            return self
        dna_seq = self.sequence.replace("U", "T")
        return AptamerSequence(
            sequence=dna_seq,
            seq_type="DNA",
            target=self.target,
            kd_nm=self.kd_nm,
            source=self.source,
            aptamer_id=self.aptamer_id,
            metadata=self.metadata,
        )

    # ------------------------------------------------------------------
    # Motif detection
    # ------------------------------------------------------------------

    def find_motif(self, motif: str) -> list[int]:
        """
        Find all start positions (0-indexed) of a motif in the sequence.

        Parameters
        ----------
        motif : str
            Nucleotide motif to search for (e.g., 'GGGG' for G-quadruplex hint).

        Returns
        -------
        list[int]
            List of start positions where motif is found.
        """
        motif = motif.upper()
        return [m.start() for m in re.finditer(f"(?={motif})", self.sequence)]

    def has_g_quadruplex_motif(self) -> bool:
        """
        Heuristically detect G-quadruplex forming motifs (G3+N1-7G3+N1-7G3+N1-7G3+).
        G-quadruplexes are associated with high-affinity aptamers.
        """
        pattern = r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
        return bool(re.search(pattern, self.sequence))

    # ------------------------------------------------------------------
    # Dunder methods
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        kd_str = f", Kd={self.kd_nm}nM" if self.kd_nm is not None else ""
        return (
            f"AptamerSequence(id='{self.aptamer_id}', "
            f"seq='{self.sequence[:20]}{'...' if self.length > 20 else ''}', "
            f"type={self.seq_type}, target='{self.target}'{kd_str})"
        )

    def __len__(self) -> int:
        return self.length

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AptamerSequence):
            return NotImplemented
        return self.sequence == other.sequence and self.seq_type == other.seq_type
