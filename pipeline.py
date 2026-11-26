import argparse
from pathlib import Path

from curation.report import write_report
from curation.stage1 import run_stage1
from curation.stage2 import run_stage2
from curation.stage3 import run_stage3
from curation.stage4 import run_stage4


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the bioactivity data curation, SAR, and reporting pipeline."
    )
    parser.add_argument(
        "input_csv",
        help="Path to the raw input CSV file (e.g., chembl_export.csv).",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default="out",
        help="Output directory for stage CSVs and report.md (default: out).",
    )
    parser.add_argument(
        "--skip-report",
        action="store_true",
        help="Run stages 1-4 only without generating report.md.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    run_stage1(args.input_csv, str(out_dir))
    run_stage2(str(out_dir / "stage1_cleaned.csv"), str(out_dir))
    run_stage3(str(out_dir / "stage2_compounds.csv"), str(out_dir))
    run_stage4(str(out_dir / "stage3_annotated.csv"), str(out_dir))

    if not args.skip_report:
        report_path = write_report(out_dir)
        print(f"Report written to {report_path}")


if __name__ == "__main__":
    main()