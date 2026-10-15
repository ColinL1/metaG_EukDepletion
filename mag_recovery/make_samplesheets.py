#!/usr/bin/env python3
"""Build nf-core/mag samplesheets from a main.nf results folder.

Reads   : <results>/trimmed_fastq/<Species>/<id>/reads/<id>_TRIM_paried_R{1,2}.fq.gz
Per-sample assemblies : <results>/assembly_based/assembly/<id>/**/<id>.contigs.fa
Co-assemblies         : <results>/coassembly_based/assembly/<group>/**/<group>.contigs.fa

group = <Species>_<method>_<buffer>, i.e. the same name main.nf gives the co-assembly.

Usage (from mag_recovery/):
    python3 make_samplesheets.py --results ../results                     # both runs
    python3 make_samplesheets.py --results ../results --run cobinning     # samplesheet.csv + samplesheet_assemblies.csv
    python3 make_samplesheets.py --results ../results --run coassembly    # samplesheet.csv + samplesheet_co_assemblies.csv
"""
import argparse
import csv
import glob
import os
import sys

READS_HEADER = ["sample", "group", "short_reads_1", "short_reads_2", "long_reads"]
ASSEMBLY_HEADER = ["id", "group", "assembler", "fasta"]


def die(msg):
    sys.exit(f"ERROR: {msg}")


def find_contigs(base, name):
    hits = glob.glob(os.path.join(base, name, "**", f"{name}.contigs.fa"), recursive=True)
    if len(hits) != 1:
        die(f"expected 1 contigs file for '{name}' under {base}, found {len(hits)}")
    return os.path.abspath(hits[0])


def collect_samples(results):
    samples = []
    for reads_dir in sorted(glob.glob(os.path.join(results, "trimmed_fastq", "*", "*", "reads"))):
        sample_dir = os.path.dirname(reads_dir)
        species = os.path.basename(os.path.dirname(sample_dir))
        sample = os.path.basename(sample_dir)
        tokens = sample.split("_")  # {species}_{replicate}_{method}_{buffer}, as in main.nf
        if len(tokens) < 4:
            die(f"sample id '{sample}' does not follow species_replicate_method_buffer")
        r1, r2 = (os.path.abspath(os.path.join(reads_dir, f"{sample}_TRIM_paried_R{i}.fq.gz")) for i in (1, 2))
        samples.append({"sample": sample, "group": f"{species}_{tokens[2]}_{tokens[3]}", "r1": r1, "r2": r2})
    if not samples:
        die(f"no trimmed reads found under {results}/trimmed_fastq/")
    return samples


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, help="main.nf --outdir")
    parser.add_argument("--out", default="samplesheets", help="output folder (default: samplesheets)")
    parser.add_argument("--run", choices=["cobinning", "coassembly"],
                        help="only write the sheets for this nf-core/mag run (default: both)")
    parser.add_argument("--no-check", action="store_true",
                        help="do not require the listed files to exist (e.g. published symlinks into another machine's work dir)")
    args = parser.parse_args()

    samples = collect_samples(args.results)
    runs = [args.run] if args.run else ["cobinning", "coassembly"]
    sheets = {"samplesheet.csv": (READS_HEADER, [[s["sample"], s["group"], s["r1"], s["r2"], ""] for s in samples])}

    if "cobinning" in runs:
        base = os.path.join(args.results, "assembly_based", "assembly")
        sheets["samplesheet_assemblies.csv"] = (ASSEMBLY_HEADER, [
            [s["sample"], s["group"], "MEGAHIT", find_contigs(base, s["sample"])] for s in samples])

    if "coassembly" in runs:
        base = os.path.join(args.results, "coassembly_based", "assembly")
        groups = sorted({s["group"] for s in samples})
        coassemblies = sorted(os.path.basename(p) for p in glob.glob(os.path.join(base, "*")) if os.path.isdir(p))
        orphans = sorted(set(coassemblies) - set(groups))
        if orphans:
            die(f"co-assemblies with no matching read group: {', '.join(orphans)}")
        sheets["samplesheet_co_assemblies.csv"] = (ASSEMBLY_HEADER, [
            [f"group-{g}", g, "MEGAHIT", find_contigs(base, g)] for g in groups])

    if not args.no_check:
        missing = [p for _, rows in sheets.values() for row in rows for p in row if p.startswith("/") and not os.path.exists(p)]
        if missing:
            die(f"{len(missing)} listed files do not exist, e.g. {missing[0]} (use --no-check to skip)")

    os.makedirs(args.out, exist_ok=True)
    for name, (header, rows) in sheets.items():
        write_csv(os.path.join(args.out, name), header, rows)


if __name__ == "__main__":
    main()
