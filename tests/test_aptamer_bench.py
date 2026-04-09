"""
Test suite for aptamer-bench.
Run with: pytest tests/ -v
"""

import pytest
from aptamer_bench.sequence import AptamerSequence
from aptamer_bench.sequence_metrics import SequenceMetrics, DiversityMetrics
from aptamer_bench.aptabase import AptaBaseLoader, BenchmarkDataset


# ── AptamerSequence tests ─────────────────────────────────────────────────────

class TestAptamerSequence:

    def test_valid_dna(self):
        apt = AptamerSequence("ACGT", seq_type="DNA")
        assert apt.sequence == "ACGT"
        assert apt.length == 4

    def test_valid_rna(self):
        apt = AptamerSequence("ACGU", seq_type="RNA")
        assert apt.seq_type == "RNA"

    def test_invalid_char_raises(self):
        with pytest.raises(ValueError):
            AptamerSequence("ACGX", seq_type="DNA")

    def test_gc_content(self):
        apt = AptamerSequence("GGCC", seq_type="DNA")
        assert apt.gc_content == 1.0

    def test_gc_content_mixed(self):
        apt = AptamerSequence("AATT", seq_type="DNA")
        assert apt.gc_content == 0.0

    def test_to_rna(self):
        apt = AptamerSequence("ACGT", seq_type="DNA")
        rna = apt.to_rna()
        assert rna.sequence == "ACGU"
        assert rna.seq_type == "RNA"

    def test_to_dna(self):
        apt = AptamerSequence("ACGU", seq_type="RNA")
        dna = apt.to_dna()
        assert dna.sequence == "ACGT"
        assert dna.seq_type == "DNA"

    def test_find_motif(self):
        apt = AptamerSequence("GGGGACGGGG", seq_type="DNA")
        positions = apt.find_motif("GGG")
        assert 0 in positions

    def test_g_quadruplex_detection(self):
        apt = AptamerSequence("GGGTGGGTGGGTGGG", seq_type="DNA")
        assert apt.has_g_quadruplex_motif() is True

    def test_lowercase_normalization(self):
        apt = AptamerSequence("acgt", seq_type="DNA")
        assert apt.sequence == "ACGT"

    def test_equality(self):
        a = AptamerSequence("ACGT", seq_type="DNA")
        b = AptamerSequence("ACGT", seq_type="DNA")
        assert a == b

    def test_nucleotide_frequencies_sum_to_one(self):
        apt = AptamerSequence("AACGT", seq_type="DNA")
        freqs = apt.nucleotide_frequencies
        assert abs(sum(freqs.values()) - 1.0) < 1e-9


# ── SequenceMetrics tests ─────────────────────────────────────────────────────

class TestSequenceMetrics:

    def setup_method(self):
        self.sm = SequenceMetrics()
        self.apt = AptamerSequence(
            "GCCTGTTGTGAGCCTCCTGTCGAA", seq_type="DNA", aptamer_id="TEST-001"
        )

    def test_compute_returns_result(self):
        result = self.sm.compute(self.apt)
        assert result.aptamer_id == "TEST-001"
        assert 0 <= result.gc_content <= 1
        assert 0 <= result.linguistic_complexity <= 1

    def test_linguistic_complexity_range(self):
        result = self.sm.compute(self.apt)
        assert 0 <= result.linguistic_complexity <= 1

    def test_batch_compute(self):
        apts = [self.apt, AptamerSequence("GGTTGGTGTGGTTGG", "DNA", aptamer_id="TEST-002")]
        results = self.sm.compute_batch(apts)
        assert len(results) == 2

    def test_shannon_entropy_uniform(self):
        # All 4 bases equally — max entropy = 2 bits
        seq = "ACGTACGTACGT"
        entropy = SequenceMetrics.shannon_entropy(seq)
        assert abs(entropy - 2.0) < 0.01

    def test_shannon_entropy_zero(self):
        # All same base — entropy = 0
        seq = "AAAAAAAAAA"
        entropy = SequenceMetrics.shannon_entropy(seq)
        assert entropy == 0.0


# ── DiversityMetrics tests ────────────────────────────────────────────────────

class TestDiversityMetrics:

    def setup_method(self):
        self.dm = DiversityMetrics()
        self.apts = [
            AptamerSequence("GCCTGTTGTGAG", "DNA"),
            AptamerSequence("ATCGATCGATCG", "DNA"),
            AptamerSequence("TTTTTTTTTTTT", "DNA"),
        ]

    def test_redundancy_all_unique(self):
        assert self.dm.redundancy_rate(self.apts) == 0.0

    def test_redundancy_all_identical(self):
        dupes = [AptamerSequence("AAAA", "DNA")] * 4
        assert self.dm.redundancy_rate(dupes) == 0.75

    def test_mean_pairwise_distance_positive(self):
        d = self.dm.mean_pairwise_distance(self.apts)
        assert d > 0

    def test_mean_pairwise_distance_single_seq(self):
        d = self.dm.mean_pairwise_distance([self.apts[0]])
        assert d == 0.0

    def test_levenshtein_identical(self):
        assert DiversityMetrics._levenshtein("ACGT", "ACGT") == 0

    def test_levenshtein_one_sub(self):
        assert DiversityMetrics._levenshtein("ACGT", "ACGG") == 1

    def test_summary_keys(self):
        summary = self.dm.summary(self.apts)
        assert "n_sequences" in summary
        assert "redundancy_rate" in summary
        assert "gc_content" in summary


# ── AptaBaseLoader tests ──────────────────────────────────────────────────────

class TestAptaBaseLoader:

    def setup_method(self):
        self.loader = AptaBaseLoader()

    def test_load_bundled(self):
        apts = self.loader.load_bundled()
        assert len(apts) > 0
        assert all(isinstance(a, AptamerSequence) for a in apts)

    def test_load_bundled_with_filter(self):
        apts = self.loader.load_bundled(target_filter="Thrombin")
        assert all("thrombin" in a.target.lower() for a in apts)

    def test_to_dataframe_columns(self):
        apts = self.loader.load_bundled()
        df = self.loader.to_dataframe(apts)
        assert "sequence" in df.columns
        assert "gc_content" in df.columns
        assert "kd_nm" in df.columns

    def test_to_fasta(self, tmp_path):
        apts = self.loader.load_bundled()
        out = tmp_path / "test.fasta"
        self.loader.to_fasta(apts, out)
        content = out.read_text()
        assert content.startswith(">")

    def test_to_json(self, tmp_path):
        import json
        apts = self.loader.load_bundled()
        out = tmp_path / "test.json"
        self.loader.to_json(apts, out)
        data = json.loads(out.read_text())
        assert isinstance(data, list)
        assert "sequence" in data[0]


# ── BenchmarkDataset tests ────────────────────────────────────────────────────

class TestBenchmarkDataset:

    def setup_method(self):
        loader = AptaBaseLoader()
        refs = loader.load_bundled(target_filter="Thrombin")
        self.bench = BenchmarkDataset(refs, name="Thrombin")
        self.candidates = [
            AptamerSequence("GGTTGGTGTGGTTGG", "DNA"),
            AptamerSequence("GCCTGTTGTGAGCCTCCTGTCGAA", "DNA"),
        ]

    def test_evaluate_returns_result(self):
        result = self.bench.evaluate(self.candidates, algorithm_name="TestAlgo")
        assert result.algorithm_name == "TestAlgo"
        assert result.n_generated == 2

    def test_compare_algorithms(self):
        outputs = {
            "AlgoA": self.candidates,
            "AlgoB": [AptamerSequence("AAAACCCCGGGGTTTT", "DNA")],
        }
        df = self.bench.compare_algorithms(outputs)
        assert len(df) == 2
        assert "algorithm" in df.columns
