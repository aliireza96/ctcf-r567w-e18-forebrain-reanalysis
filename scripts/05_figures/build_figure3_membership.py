#!/usr/bin/env python3
"""Freeze the direct GO and curated memberships used by Figure 3B-D.

Usage:
  python build_figure3_membership.py ORG_MM_EG_SQLITE [OUTPUT_CSV]

The SQLite file must be from org.Mm.eg.db 3.20.0 (Entrez/GO source date
2024-Sep20).  The generated supplementary CSV is self-contained; downstream
analysis reads that CSV and does not require R or an annotation database.
"""

import argparse
import csv
import sqlite3
from collections import OrderedDict
from pathlib import Path


SETS = OrderedDict([
    ("Axon guidance", ("B", "Guidance", "GO:0007411", "direct GO", "0.87", 191)),
    ("Semaphorin-plexin signaling", ("B", "Guidance", "GO:0071526", "direct GO", "1.0", 43)),
    ("Cell adhesion", ("B", "Adhesion", "GO:0007155", "direct GO", "0.58", 444)),
    ("Synaptic transmission", ("C", "Neurotransmission", "GO:0007268", "direct GO", "0.59", 161)),
    ("Ionotropic glutamate receptor sig", ("C", "Neurotransmission", "GO:0035235", "direct GO", "1.0", 22)),
    ("Neuronal action potential", ("C", "Excitability", "GO:0019228", "direct GO", "0.62", 53)),
    ("Cell cycle (S and G2/M)", ("D", "Proliferation", "curated (Tirosh 2016)", "curated", "1.0", 94)),
    ("DNA repair", ("D", "DNA repair", "GO:0006281", "direct GO", "0.62", 253)),
    ("Response to oxidative stress", ("D", "Stress", "GO:0006979", "direct GO", "0.87", 144)),
    ("Apoptotic process", ("D", "Cell death", "GO:0006915", "direct GO", "0.72", 533)),
])

CELL_CYCLE = """Anln Anp32e Atad2 Aurka Aurkb Birc5 Blm Brip1 Bub1 Casp8ap2 Cbx5
Ccnb2 Ccne2 Cdc20 Cdc25c Cdc45 Cdc6 Cdca2 Cdca3 Cdca7 Cdca8 Cdk1 Cenpa Cenpe
Cenpf Chaf1b Ckap2 Ckap2l Ckap5 Cks1b Cks2 Clspn Ctcf Dlgap5 Dscc1 Dtl E2f8
Ect2 Exo1 Fen1 G2e3 Gas2l3 Gins2 Gmnn Gtse1 Hells Hjurp Hmgb2 Hmmr Kif11
Kif20b Kif23 Kif2c Lbr Mcm2 Mcm4 Mcm5 Mcm6 Mki67 Msh2 Nasp Ncapd2 Ndc80
Nek2 Nuf2 Nusap1 Pcna Pola1 Pold3 Prim1 Psrc1 Rad51 Rad51ap1 Rangap1 Rfc2
Rpa2 Rrm1 Rrm2 Slbp Smc4 Tacc3 Tipin Tmpo Top2a Tpx2 Ttk Tubb4b Tyms Ube2c
Ubr7 Uhrf1 Ung Usp1 Wdr76""".split()


def metadata(connection):
    return dict(connection.execute("SELECT name, value FROM metadata"))


def direct_genes(connection, go_id):
    query = """
        SELECT DISTINCT gene_info.symbol
        FROM go_bp JOIN gene_info USING (_id)
        WHERE go_bp.go_id = ?
        ORDER BY gene_info.symbol
    """
    return [row[0] for row in connection.execute(query, (go_id,))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sqlite", type=Path)
    parser.add_argument("output", nargs="?", type=Path,
                        default=Path(__file__).resolve().parents[2] / "tables" /
                                "Table_S11_gene_set_membership.csv")
    args = parser.parse_args()
    connection = sqlite3.connect(args.sqlite.resolve())
    db_meta = metadata(connection)
    if db_meta.get("ORGANISM") != "Mus musculus" or db_meta.get("GOEGSOURCEDATE") != "2024-Sep20":
        raise ValueError("Expected org.Mm.eg.db 3.20.0 / GOEGSOURCEDATE 2024-Sep20")

    memberships = OrderedDict()
    for name, (_panel, _category, go_id, annotation, _fraction, expected) in SETS.items():
        genes = sorted(CELL_CYCLE) if annotation == "curated" else direct_genes(connection, go_id)
        if len(genes) != expected:
            raise ValueError(f"{name}: expected {expected} genes, found {len(genes)}")
        memberships[name] = genes

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        handle.write("# Authoritative membership for the 10 gene sets displayed in Figure 3B-D.\n")
        handle.write("# Direct rows use org.Mm.eg.db 3.20.0, keytype='GO' (not GOALL),\n")
        handle.write("# Entrez/GO source date 2024-Sep20. The cell-cycle row is the recovered\n")
        handle.write("# 94-gene mouse-symbol Tirosh 2016 S+G2/M set used in the analysis.\n")
        handle.write("# One row per unique set-member gene; this table is the frozen input to\n")
        handle.write("# compute_figure3_tables.py and removes annotation-release ambiguity.\n#\n")
        fields = ["panel", "gene_set", "category", "go_id", "annotation",
                  "canonical_fraction", "source_release", "gene"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name, genes in memberships.items():
            panel, category, go_id, annotation, fraction, _expected = SETS[name]
            source = ("Tirosh 2016 S+G2/M; recovered mouse-symbol set"
                      if annotation == "curated" else
                      "org.Mm.eg.db 3.20.0; GOEGSOURCEDATE 2024-Sep20")
            for gene in genes:
                writer.writerow({"panel": panel, "gene_set": name, "category": category,
                                 "go_id": go_id, "annotation": annotation,
                                 "canonical_fraction": fraction,
                                 "source_release": source, "gene": gene})
    print(f"wrote {args.output.resolve()} ({sum(map(len, memberships.values()))} memberships)")


if __name__ == "__main__":
    main()
