"""
aptamer-bench Streamlit Dashboard

Run with: streamlit run aptamer_bench/visualization/dashboard.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

try:
    import logomaker
    HAS_LOGOMAKER = True
except ImportError:
    HAS_LOGOMAKER = False

from aptamer_bench.aptabase import AptaBaseLoader
from benchmark import BenchmarkDataset
from aptamer_bench.sequence_metrics import SequenceMetrics
from diversity_metrics import DiversityMetrics
from aptamer_bench.sequence import AptamerSequence

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="aptamer-bench",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;700&family=Inter:wght@400;600&display=swap');
  html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
  .stMetric label { font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.05em; }
  .block-container { padding-top: 1.5rem; }
  h1, h2, h3 { font-family: 'JetBrains Mono', monospace; }
</style>
""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🧬 aptamer-bench")
    st.caption("Benchmarking toolkit for in-silico aptamer evaluation")
    st.divider()

    mode = st.radio(
        "Mode",
        ["📊 Dataset Explorer", "🔬 Sequence Analyzer", "⚖️ Algorithm Comparison"],
    )
    st.divider()
    st.markdown("**Compatible with** [pyaptamer](https://github.com/gc-os-ai/pyaptamer)")

# ── Data loading ──────────────────────────────────────────────────────────────
@st.cache_data
def load_bundled_data():
    loader = AptaBaseLoader()
    aptamers = loader.load_bundled()
    df = loader.to_dataframe(aptamers)
    return aptamers, df

aptamers, df = load_bundled_data()

# ═════════════════════════════════════════════════════════════════════════════
# MODE 1 — Dataset Explorer
# ═════════════════════════════════════════════════════════════════════════════
if mode == "📊 Dataset Explorer":
    st.title("Dataset Explorer")
    st.markdown("Explore the bundled curated aptamer dataset. Export to CSV, JSON, or FASTA.")

    # Filters
    col1, col2 = st.columns([2, 1])
    with col1:
        target_options = ["All"] + sorted(df["target"].unique().tolist())
        selected_target = st.selectbox("Filter by target", target_options)
    with col2:
        seq_type_filter = st.selectbox("Sequence type", ["All", "DNA", "RNA"])

    filtered = df.copy()
    if selected_target != "All":
        filtered = filtered[filtered["target"] == selected_target]
    if seq_type_filter != "All":
        filtered = filtered[filtered["seq_type"] == seq_type_filter]

    # KPIs
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Sequences", len(filtered))
    c2.metric("Unique targets", filtered["target"].nunique())
    c3.metric(
        "Median Kd (nM)",
        f"{filtered['kd_nm'].median():.2f}" if filtered["kd_nm"].notna().any() else "N/A",
    )
    c4.metric("Mean GC content", f"{filtered['gc_content'].mean():.2%}")

    st.dataframe(filtered, use_container_width=True, hide_index=True)

    # GC distribution plot
    st.subheader("GC Content Distribution")
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.hist(filtered["gc_content"], bins=10, color="#4F8EF7", edgecolor="white", alpha=0.85)
    ax.set_xlabel("GC Content")
    ax.set_ylabel("Count")
    ax.spines[["top", "right"]].set_visible(False)
    st.pyplot(fig, use_container_width=True)

    # Export
    st.subheader("Export")
    col_e1, col_e2 = st.columns(2)
    with col_e1:
        csv = filtered.to_csv(index=False)
        st.download_button("⬇ Download CSV", csv, "aptamers.csv", "text/csv")
    with col_e2:
        fasta_lines = []
        for _, row in filtered.iterrows():
            fasta_lines.append(
                f">{row['aptamer_id']} target={row['target']} kd_nm={row['kd_nm']}\n{row['sequence']}"
            )
        st.download_button(
            "⬇ Download FASTA",
            "\n".join(fasta_lines),
            "aptamers.fasta",
            "text/plain",
        )

# ═════════════════════════════════════════════════════════════════════════════
# MODE 2 — Sequence Analyzer
# ═════════════════════════════════════════════════════════════════════════════
elif mode == "🔬 Sequence Analyzer":
    st.title("Sequence Analyzer")
    st.markdown("Paste one or more aptamer sequences to compute quality metrics.")

    default_seqs = "GGTTGGTGTGGTTGG\nGCCTGTTGTGAGCCTCCTGTCGAA\nATCCGATCGATCGAT"
    raw_input = st.text_area(
        "Enter sequences (one per line)",
        value=default_seqs,
        height=150,
    )
    seq_type_input = st.radio("Sequence type", ["DNA", "RNA"], horizontal=True)

    if st.button("Analyze", type="primary"):
        lines = [l.strip().upper() for l in raw_input.strip().splitlines() if l.strip()]
        parsed = []
        errors = []
        for i, seq in enumerate(lines):
            try:
                parsed.append(
                    AptamerSequence(seq, seq_type=seq_type_input, aptamer_id=f"INPUT-{i+1:02d}")
                )
            except ValueError as e:
                errors.append(f"Line {i+1}: {e}")

        if errors:
            for err in errors:
                st.error(err)

        if parsed:
            sm = SequenceMetrics()
            dm = DiversityMetrics()
            results = sm.compute_batch(parsed)
            metrics_df = pd.DataFrame([r.to_dict() for r in results])

            st.subheader("Per-Sequence Metrics")
            st.dataframe(metrics_df, use_container_width=True, hide_index=True)

            st.subheader("Population Diversity")
            div_summary = dm.summary(parsed)
            c1, c2, c3 = st.columns(3)
            c1.metric("Unique sequences", div_summary["n_unique"])
            c2.metric("Redundancy rate", f"{div_summary['redundancy_rate']:.2%}")
            c3.metric(
                "Mean pairwise distance",
                f"{div_summary['mean_pairwise_distance_normalized']:.3f}",
            )

            if HAS_LOGOMAKER and len(parsed) >= 2:
                # Sequence logo — only if all same length
                lengths = {len(a.sequence) for a in parsed}
                if len(lengths) == 1:
                    st.subheader("Sequence Logo")
                    bases = "ACGT" if seq_type_input == "DNA" else "ACGU"
                    L = list(lengths)[0]
                    counts = {b: [0] * L for b in bases}
                    for apt in parsed:
                        for j, nuc in enumerate(apt.sequence):
                            if nuc in counts:
                                counts[nuc][j] += 1
                    logo_df = pd.DataFrame(counts)
                    logo_df = logo_df.div(logo_df.sum(axis=1), axis=0)
                    fig, ax = plt.subplots(figsize=(max(6, L * 0.5), 2.5))
                    logomaker.Logo(logo_df, ax=ax, color_scheme="classic")
                    ax.set_xlabel("Position")
                    st.pyplot(fig, use_container_width=True)

# ═════════════════════════════════════════════════════════════════════════════
# MODE 3 — Algorithm Comparison
# ═════════════════════════════════════════════════════════════════════════════
elif mode == "⚖️ Algorithm Comparison":
    st.title("Algorithm Comparison")
    st.markdown(
        "Simulate comparing multiple in-silico aptamer generation algorithms "
        "against a reference dataset. In real use, paste outputs from pyaptamer algorithms."
    )

    st.info(
        "This demo uses synthetic data to illustrate the comparison interface. "
        "In practice, connect algorithm outputs from pyaptamer (AptaNet, AptaTrans, MAWS).",
        icon="ℹ️",
    )

    # Generate synthetic demo data
    @st.cache_data
    def make_synthetic_outputs():
        rng = np.random.default_rng(42)
        bases_dna = list("ACGT")

        def random_seqs(n, length_range, gc_bias=0.5):
            seqs = []
            for _ in range(n):
                L = rng.integers(*length_range)
                seq = "".join(
                    rng.choice(
                        bases_dna,
                        p=[
                            (1 - gc_bias) / 2,
                            gc_bias / 2,
                            gc_bias / 2,
                            (1 - gc_bias) / 2,
                        ],
                    )
                    for _ in range(L)
                )
                seqs.append(AptamerSequence(seq, "DNA", kd_nm=float(rng.uniform(0.1, 100))))
            return seqs

        return {
            "AptaNet": random_seqs(20, (18, 35), gc_bias=0.55),
            "AptaTrans": random_seqs(20, (20, 30), gc_bias=0.45),
            "MAWS": random_seqs(20, (15, 40), gc_bias=0.50),
        }

    outputs = make_synthetic_outputs()
    loader = AptaBaseLoader()
    refs = loader.load_bundled(target_filter="Thrombin")
    bench = BenchmarkDataset(refs, name="Thrombin-Demo")

    comparison_df = bench.compare_algorithms(outputs)
    st.subheader("Summary Comparison Table")
    st.dataframe(comparison_df, use_container_width=True, hide_index=True)

    # Bar charts
    st.subheader("Redundancy Rate by Algorithm")
    fig, axes = st.columns(2)

    with axes[0]:
        fig1, ax1 = plt.subplots(figsize=(5, 3))
        colors = ["#4F8EF7", "#F7844F", "#4FF7A0"]
        ax1.bar(comparison_df["algorithm"], comparison_df["redundancy_rate"], color=colors)
        ax1.set_ylabel("Redundancy Rate")
        ax1.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
        ax1.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig1, use_container_width=True)

    with axes[1]:
        fig2, ax2 = plt.subplots(figsize=(5, 3))
        ax2.bar(comparison_df["algorithm"], comparison_df["mean_pairwise_distance"], color=colors)
        ax2.set_ylabel("Mean Pairwise Distance (norm.)")
        ax2.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig2, use_container_width=True)

    st.download_button(
        "⬇ Download Comparison CSV",
        comparison_df.to_csv(index=False),
        "algorithm_comparison.csv",
        "text/csv",
    )
