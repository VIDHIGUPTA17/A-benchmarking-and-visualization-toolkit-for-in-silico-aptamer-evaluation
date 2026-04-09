# aptamer-bench 🧬

**A benchmarking and visualization toolkit for in-silico aptamer evaluation.**

Designed to complement [`pyaptamer`](https://github.com/gc-os-ai/pyaptamer) by providing
standardized benchmark datasets, evaluation metrics, and visualization tools — directly
addressing the gap stated in the pyaptamer roadmap:

> *"There are no publicly available, standardized benchmark datasets or challenges
> to test the performance of in-silico aptamer design procedures."*

---

## Features

| Module | Description |
|---|---|
| `AptamerSequence` | Core sequence object (DNA/RNA) with GC content, G-quadruplex detection, motif search |
| `SequenceMetrics` | Linguistic complexity, Shannon entropy, structural hints |
| `DiversityMetrics` | Levenshtein pairwise distance, redundancy rate, population GC/length stats |
| `AptaBaseLoader` | Load curated aptamers from bundled dataset or CSV; export to FASTA/JSON/CSV |
| `BenchmarkDataset` | Compare multiple algorithms side-by-side with standardized metrics |
| Streamlit Dashboard | Interactive visualization of sequences, metrics, and algorithm comparison |

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/aptamer-bench.git
cd aptamer-bench
pip install -e ".[dev]"
```

---

## Quick Start

```python
from aptamer_bench import AptamerSequence, SequenceMetrics, AptaBaseLoader, BenchmarkDataset

# 1. Create a sequence
apt = AptamerSequence("GGTTGGTGTGGTTGG", seq_type="DNA", target="Thrombin", kd_nm=26.0)
print(apt.gc_content)          # 0.533
print(apt.has_g_quadruplex_motif())  # True

# 2. Load curated dataset
loader = AptaBaseLoader()
aptamers = loader.load_bundled(target_filter="Thrombin")
loader.to_fasta(aptamers, "thrombin.fasta")   # pyaptamer-compatible

# 3. Benchmark algorithm outputs
bench = BenchmarkDataset(aptamers, name="Thrombin")
my_candidates = [AptamerSequence("GCCTGTTGTGAG", "DNA")]
result = bench.evaluate(my_candidates, algorithm_name="MyAlgo")
print(result)
```

Run the full example:
```bash
python quickstart.py
```

---

## Streamlit Dashboard

```bash
streamlit run dashboard.py
```

Features:
- **Dataset Explorer** — browse, filter, and export curated aptamers
- **Sequence Analyzer** — paste sequences and compute metrics instantly
- **Algorithm Comparison** — compare pyaptamer algorithms (AptaNet, AptaTrans, MAWS)

---

## Running Tests

```bash
pytest tests/ -v --cov=aptamer_bench
```

---

## Roadmap

- [ ] Live integration with AptaBase / aptamer.atdbio.com API
- [ ] Secondary structure prediction (ViennaRNA wrapper)
- [ ] Binding affinity regression baseline model
- [ ] pyaptamer algorithm output parsers (AptaNet, AptaTrans, MAWS)
- [ ] Standardized challenge dataset publication

---

## Relation to pyaptamer

`aptamer-bench` is designed as a benchmarking companion to
[pyaptamer](https://github.com/gc-os-ai/pyaptamer). It exports data in formats
directly consumable by pyaptamer algorithms, and accepts their outputs for evaluation.

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Contributing

Issues and PRs welcome! This project was created to support the
[ESoC 2026 pyaptamer project](https://github.com/gc-os-ai/pyaptamer).
