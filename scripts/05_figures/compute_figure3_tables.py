#!/usr/bin/env python3
"""Recompute the authoritative Figure 3B-D gene-set battery.

The analysis uses the frozen membership in Table_S11_gene_set_membership.csv,
4,000 size-matched permutation draws within each cell-class DE universe, and a
15-tested-gene floor.  It writes the complete 10-set x 9-class display to
Table_S3_gene_set_battery.csv.  Cells below the floor remain in the rectangular
display but are explicitly NA and are excluded from Benjamini-Hochberg correction.

Positional arguments are optional and, in order, are DE_DIR, MEMBERSHIP_CSV and
OUTPUT_CSV.  Defaults are resolved from the script location, not the working
directory.  Requires numpy.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import OrderedDict
from pathlib import Path

import numpy as np


N_DRAWS = 4000
RNG_SEED = 3
MIN_TESTED_GENES = 15
PREFIX = "05b_cortstr_plus_hypoglut_"

CLASSES = [
    ("Cycling RG", "DE_celltype_broad_Cycling_RG", "Dorsal"),
    ("Immature astrocytes", "DE_celltype_broad_Immature_Astrocytes", "Dorsal"),
    ("Deep-layer EN", "DE_celltype_broad_Deep_layer_EN", "Dorsal"),
    ("Upper-layer EN", "DE_celltype_broad_Upper_layer_EN", "Dorsal"),
    ("SPN-D1", "DE_celltype_fine_Striatal_SPN_D1_striosome_like_", "LGE-derived"),
    ("SPN-D2", "DE_celltype_fine_Striatal_SPN_D2_indirect_pathway_", "LGE-derived"),
    ("LGE-IN precursors", "DE_celltype_broad_LGE_IN_prec", "LGE-derived"),
    ("MGE interneurons", "DE_celltype_broad_MGE_IN", "Cortical IN"),
    ("CGE interneurons", "DE_celltype_broad_Migrating_CGE_derived_IN", "Cortical IN"),
]

EXPECTED_SET_SIZES = OrderedDict([
    ("Axon guidance", 191),
    ("Semaphorin-plexin signaling", 43),
    ("Cell adhesion", 444),
    ("Synaptic transmission", 161),
    ("Ionotropic glutamate receptor sig", 22),
    ("Neuronal action potential", 53),
    ("Cell cycle (S and G2/M)", 94),
    ("DNA repair", 253),
    ("Response to oxidative stress", 144),
    ("Apoptotic process", 533),
])


def comment_aware_rows(path: Path):
    """Read a CSV while ignoring blank and comment-prefixed records."""
    with path.open(newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


def locate_defaults(script: Path):
    table_candidates = [script.parent.parent / "tables", script.parents[2] / "tables"]
    tables = next((p for p in table_candidates if (p / "Table_S11_gene_set_membership.csv").exists()),
                  table_candidates[0])
    project = next((p for p in script.parents
                    if (p / "E18p5_clean" / "results" / "tables").is_dir()), None)
    if project is None:
        raise FileNotFoundError("Could not locate E18p5_clean/results/tables from the script path")
    return project / "E18p5_clean" / "results" / "tables", tables


def load_membership(path: Path):
    rows = comment_aware_rows(path)
    required = {"panel", "gene_set", "category", "go_id", "annotation",
                "canonical_fraction", "gene"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} lacks required membership columns: {sorted(required)}")

    sets = OrderedDict()
    for row in rows:
        name = row["gene_set"]
        meta = tuple(row[k] for k in
                     ("panel", "category", "go_id", "annotation", "canonical_fraction"))
        if name not in sets:
            sets[name] = {"meta": meta, "genes": []}
        elif sets[name]["meta"] != meta:
            raise ValueError(f"Inconsistent metadata for {name}")
        sets[name]["genes"].append(row["gene"])

    if list(sets) != list(EXPECTED_SET_SIZES):
        raise ValueError("Membership set names/order do not match the authoritative 10-set battery")
    for name, expected in EXPECTED_SET_SIZES.items():
        genes = sets[name]["genes"]
        if len(genes) != len(set(genes)):
            raise ValueError(f"Duplicate genes in {name}")
        if len(genes) != expected:
            raise ValueError(f"{name}: expected {expected} genes, found {len(genes)}")
    return sets


def load_de(path: Path):
    rows = comment_aware_rows(path)
    out = OrderedDict()
    for row in rows:
        try:
            out[row["gene"]] = float(row["avg_log2FC"])
        except (KeyError, TypeError, ValueError):
            continue
    if not out:
        raise ValueError(f"No gene/avg_log2FC values read from {path}")
    return out


def test_rng(master_seed: int, gene_set: str, cell_type: str):
    digest = hashlib.sha256(f"{master_seed}|{gene_set}|{cell_type}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def permutation_test(values, member_genes, rng, draws):
    members = [g for g in member_genes if g in values]
    k = len(members)
    all_values = np.asarray(list(values.values()), dtype=float)
    observed = float(np.mean([values[g] for g in members]))
    null = np.empty(draws, dtype=float)
    for index in range(draws):
        null[index] = rng.choice(all_values, size=k, replace=False).mean()
    centre = float(null.mean())
    spread = float(null.std(ddof=0))
    z_score = (observed - centre) / spread if spread else 0.0
    p_value = (int(np.count_nonzero(np.abs(null - centre) >= abs(observed - centre))) + 1) / (draws + 1)
    return observed, z_score, p_value, k


def bh_adjust(p_values):
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out.tolist()


def binomial_tail(successes, trials, probability=0.05):
    return sum(math.comb(trials, k) * probability**k * (1 - probability)**(trials - k)
               for k in range(successes, trials + 1))


def fmt(value, digits):
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def compute(de_dir: Path, membership_path: Path, draws: int, seed: int, floor: int):
    sets = load_membership(membership_path)
    de = {}
    for cell_type, stem, _group in CLASSES:
        de[cell_type] = load_de(de_dir / f"{PREFIX}{stem}_mut_vs_wt.csv")

    rows, tested_indices, p_values = [], [], []
    for gene_set, payload in sets.items():
        panel, category, go_id, annotation, canonical_fraction = payload["meta"]
        genes = payload["genes"]
        for cell_type, _stem, group in CLASSES:
            n_tested = sum(g in de[cell_type] for g in genes)
            row = OrderedDict([
                ("panel", panel), ("gene_set", gene_set), ("category", category),
                ("go_id", go_id), ("annotation", annotation),
                ("canonical_fraction", canonical_fraction),
                ("cell_type", cell_type), ("group", group),
                ("n_genes_in_set", str(len(genes))), ("n_genes_tested", str(n_tested)),
            ])
            if n_tested < floor:
                row.update(test_status=f"NA_below_{floor}_genes", mean_log2FC="NA",
                           z_vs_random="NA", perm_p="NA", significant="NA")
            else:
                obs, z_score, p_value, _ = permutation_test(
                    de[cell_type], genes, test_rng(seed, gene_set, cell_type), draws)
                row.update(test_status="tested", mean_log2FC=fmt(obs, 4),
                           z_vs_random=fmt(z_score, 3), perm_p=fmt(p_value, 4),
                           significant=str(p_value < 0.05))
                tested_indices.append(len(rows))
                p_values.append(p_value)
            rows.append(row)

    q_values = bh_adjust(p_values)
    for row in rows:
        row["significant_FDR_q05"] = "NA"
        row["q_value_BH"] = "NA"
    for index, q_value in zip(tested_indices, q_values):
        rows[index]["significant_FDR_q05"] = str(q_value < 0.05)
        rows[index]["q_value_BH"] = fmt(q_value, 4)
    return rows, len(p_values)


def write_table(path: Path, rows, tested_count: int, draws: int, seed: int, floor: int):
    nominal = sum(row["significant"] == "True" for row in rows)
    fdr = sum(row["significant_FDR_q05"] == "True" for row in rows)
    numeric_q = [float(row["q_value_BH"]) for row in rows if row["q_value_BH"] != "NA"]
    tail = binomial_tail(nominal, tested_count)
    comments = [
        "# Authoritative source for Figure 3B-D: 10 displayed gene sets x 9 cell classes.",
        "# Membership is frozen in Table_S11_gene_set_membership.csv: direct GO annotations",
        "# from org.Mm.eg.db 3.20.0 (GO/Entrez source date 2024-Sep20), plus the curated",
        "# 94-gene Tirosh S+G2/M cell-cycle set. GOALL/inherited annotations are not used.",
        f"# Each test compares its mean log2FC with {draws:,} size-matched random gene sets",
        f"# sampled from the same cell-class DE universe (deterministic keyed seed {seed}).",
        f"# A hard floor of {floor} tested genes is applied per set/class cell.",
        "# The cell-cycle set is retained as a targeted biological-control set even though it",
        "# is not testable in every class; its seven below-floor cells are explicitly NA.",
        f"# Therefore the rectangular display has 90 cells but only {tested_count} inferential tests.",
        f"# Benjamini-Hochberg q values are calculated across those {tested_count} non-NA tests only.",
        f"# Nominal P<0.05 cells: {nominal}/{tested_count}; exact Binomial({tested_count},0.05) upper-tail P={tail:.3g}.",
        f"# FDR q<0.05 cells: {fdr}/{tested_count}; lowest q={min(numeric_q):.4f}.",
        "# Panels are descriptive and carry no significance symbols; per-cell P and q values",
        "# are retained for transparency. Removed DNA-replication and cytoplasmic-translation",
        "# sets remain documented in Table_S6_gene_set_audit.csv and are not part of this table.",
        "#",
    ]
    fieldnames = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write("\n".join(comments) + "\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return nominal, fdr, min(numeric_q), tail


def main():
    script = Path(__file__).resolve()
    default_de, default_tables = locate_defaults(script)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("de_dir", nargs="?", type=Path, default=default_de)
    parser.add_argument("membership_csv", nargs="?", type=Path,
                        default=default_tables / "Table_S11_gene_set_membership.csv")
    parser.add_argument("output_csv", nargs="?", type=Path,
                        default=default_tables / "Table_S3_gene_set_battery.csv")
    parser.add_argument("--draws", type=int, default=N_DRAWS)
    parser.add_argument("--seed", type=int, default=RNG_SEED)
    parser.add_argument("--min-tested-genes", type=int, default=MIN_TESTED_GENES)
    args = parser.parse_args()
    rows, tested = compute(args.de_dir.resolve(), args.membership_csv.resolve(),
                           args.draws, args.seed, args.min_tested_genes)
    nominal, fdr, lowest, tail = write_table(args.output_csv.resolve(), rows, tested,
                                              args.draws, args.seed, args.min_tested_genes)
    print(f"wrote {args.output_csv.resolve()}")
    print(f"display cells=90; inferential tests={tested}; below-floor NA={90-tested}")
    print(f"nominal P<0.05={nominal}; FDR q<0.05={fdr}; lowest q={lowest:.4f}; binomial P={tail:.3g}")


if __name__ == "__main__":
    main()
