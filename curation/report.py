from pathlib import Path
import pandas as pd

from curation.config import (
    DELTA_THRESHOLD,
    FINGERPRINT_BITS,
    FINGERPRINT_RADIUS,
    HIGH_VARIANCE_RANGE,
    SIM_THRESHOLD,
)

CLIFF_COLS = ["chembl_id_a", "chembl_id_b", "pIC50_a", "pIC50_b",
              "similarity", "delta_pIC50", "same_scaffold"]


def _counts(series: pd.Series, n: int | None = None) -> dict:
    c = series.dropna().value_counts()
    return {k: int(v) for k, v in (c.head(n) if n else c).items()}


def collect_stats(out_dir: str | Path) -> dict:
    p = Path(out_dir)
    rej = pd.read_csv(p / "rejected.csv")
    s1 = pd.read_csv(p / "stage1_cleaned.csv")
    s2 = pd.read_csv(p / "stage2_compounds.csv")
    s3 = pd.read_csv(p / "stage3_annotated.csv")
    s4 = pd.read_csv(p / "stage4_annotated.csv")
    cl = pd.read_csv(p / "activity_cliffs.csv")

    exact = s2[s2["activity_type"] == "exact"]
    cens = s2[s2["activity_type"] == "censored_only"]
    n_exact_rows = int((~s1["censored"]).sum())
    scaf = s4["scaffold"].value_counts()

    return {
        "targets": sorted(s1["Target Name"].dropna().unique()),
        "n_input": len(s1) + len(rej),
        "n_cleaned": len(s1),
        "n_rejected": len(rej),
        "reject_reasons": _counts(rej["reject_reason"]),
        "n_had_salt": int(s1["had_salt"].sum()),
        "n_chembl_ids": int(s1["Molecule ChEMBL ID"].nunique()),
        "n_parent_structures": int(s1["parent_inchikey"].nunique()),
        "n_compounds": len(s2),
        "n_exact": len(exact),
        "n_censored_only": len(cens),
        "n_high_variance": int(s2["high_variance"].sum()),
        "n_exact_rows": n_exact_rows,
        "n_duplicates_merged": n_exact_rows - len(exact),
        "n_single_measurement": int((exact["n_measurements"] == 1).sum()),
        "censored_relations": _counts(cens["pIC50_relation"]),
        "lipinski_counts": {int(k): int(v) for k, v in
                            s3["lipinski_violations"].value_counts().sort_index().items()},
        "n_pains": int(s3["pains_flag"].sum()),
        "top_pains": _counts(s3.loc[s3["pains_flag"], "pains_name"], 3),
        "mw_mean": float(s3["MW"].mean()),
        "clogp_mean": float(s3["cLogP"].mean()),
        "n_scaffolds": len(scaf),
        "n_singleton_scaffolds": int((scaf == 1).sum()),
        "top_scaffolds": _counts(s4["scaffold"], 5),
        "n_cliffs": len(cl),
        "n_cliff_compounds": int(pd.concat([cl["chembl_id_a"], cl["chembl_id_b"]]).nunique()),
        "n_cliffs_same_scaffold": int(cl["same_scaffold"].sum()),
        "n_cliffs_same_connectivity": int(cl["same_connectivity"].sum()),
        "top_cliffs": cl.head(5)[CLIFF_COLS].to_dict("records"),
    }


def _table(head: list[str], rows: list[list]) -> list[str]:
    rows = rows or [["none"] + ["-"] * (len(head) - 1)]
    return (["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
            + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows])


def render_report(s: dict) -> str:
    L = ["# Bioactivity Curation & SAR Report"]

    def section(title: str, *body: str) -> None:
        L.extend(["", f"## {title}", "", *body])

    targets = ", ".join(s["targets"])
    warn = (["", "**Warning: more than one target in the input; potencies were merged across targets.**"]
            if len(s["targets"]) > 1 else [])

    section("1. Summary",
        f"- Target: {targets} (IC50, binding assays).",
        f"- {s['n_input']:,} input rows: {s['n_cleaned']:,} cleaned, {s['n_rejected']:,} rejected.",
        f"- {s['n_compounds']:,} compounds: {s['n_exact']:,} with exact potency, "
        f"{s['n_censored_only']:,} censored-only.",
        *warn)

    section("2. Rejected rows",
        *_table(["Reason", "Count"], [[k, f"{v:,}"] for k, v in s["reject_reasons"].items()]))

    section("3. Standardization",
        f"- {s['n_had_salt']:,} rows contained salts or mixtures (largest fragment kept).",
        f"- {s['n_chembl_ids']:,} ChEMBL IDs map to {s['n_parent_structures']:,} parent structures.")

    section("4. Activity data",
        f"- {s['n_exact_rows']:,} exact measurements merged into {s['n_exact']:,} compounds "
        f"({s['n_duplicates_merged']:,} duplicate measurements merged).",
        f"- {s['n_high_variance']:,} compounds flagged: measurements span more than "
        f"{HIGH_VARIANCE_RANGE} pIC50.", "",
        *_table(["Censored-only relation (pIC50)", "Count"],
                [[k, f"{v:,}"] for k, v in s["censored_relations"].items()]))

    section("5. Properties and flags",
        f"- Mean MW {s['mw_mean']:.1f}, mean cLogP {s['clogp_mean']:.2f}.",
        f"- {s['n_pains']:,} compounds match a PAINS pattern.", "",
        *_table(["Lipinski violations", "Compounds"],
                [[k, f"{v:,}"] for k, v in s["lipinski_counts"].items()]), "",
        *_table(["Top PAINS patterns", "Compounds"],
                [[k, f"{v:,}"] for k, v in s["top_pains"].items()]))

    section("6. SAR and activity cliffs",
        f"- {s['n_scaffolds']:,} Murcko scaffolds ({s['n_singleton_scaffolds']:,} singletons).",
        f"- {s['n_cliffs']:,} cliff pairs involving {s['n_cliff_compounds']:,} compounds "
        f"({s['n_cliffs_same_scaffold']:,} within one scaffold, "
        f"{s['n_cliffs_same_connectivity']:,} stereoisomer pairs).", "",
        *_table(["Top scaffolds", "Compounds"],
                [[f"`{k}`", f"{v:,}"] for k, v in s["top_scaffolds"].items()]), "",
        *_table(["A", "B", "pIC50 A", "pIC50 B", "Similarity", "Delta", "Same scaffold"],
                [[c["chembl_id_a"], c["chembl_id_b"], f"{c['pIC50_a']:.2f}", f"{c['pIC50_b']:.2f}",
                  f"{c['similarity']:.2f}", f"{c['delta_pIC50']:.2f}", c["same_scaffold"]]
                 for c in s["top_cliffs"]]))

    section("7. Parameters",
        f"- High-variance range: more than {HIGH_VARIANCE_RANGE} pIC50",
        f"- Fingerprint: Morgan, radius {FINGERPRINT_RADIUS}, {FINGERPRINT_BITS} bits",
        f"- Activity cliff: Tanimoto at least {SIM_THRESHOLD} and delta pIC50 at least {DELTA_THRESHOLD}")

    top = next(iter(s["top_pains"].items()), None)
    pains_note = (f"- {top[0]} accounts for {top[1]:,} of {s['n_pains']:,} PAINS flags; "
                  "flags are alerts, not proof." if top else "- PAINS flags are alerts, not proof.")
    section("8. Limitations",
        f"- {s['n_single_measurement']:,} compounds have one measurement, so the variance check cannot flag them.",
        f"- {s['n_censored_only']:,} censored-only compounds are excluded from variance and cliff analysis.",
        f"- {s['n_high_variance']:,} high-variance compounds are also excluded from cliff analysis.",
        "- Potencies are medians across assays, so some cliffs may reflect assay differences, not chemistry.",
        "- Morgan fingerprints ignore stereochemistry, so stereoisomers can look identical.",
        pains_note,
        f"- {s['reject_reasons'].get('non_molar_units', 0):,} rows reported in mass units (ug/mL) were rejected, not converted.")

    return "\n".join(L) + "\n"


def write_report(out_dir: str | Path) -> Path:
    path = Path(out_dir) / "report.md"
    path.write_text(render_report(collect_stats(out_dir)), encoding="utf-8")
    return path


if __name__ == "__main__":
    print(write_report("out"))