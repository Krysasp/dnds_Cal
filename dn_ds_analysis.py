#!/usr/bin/env python3
import argparse
import os
import subprocess
import itertools
import tempfile
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import pycountry
import pycountry_convert as pc
from thefuzz import process
import logging
import re
import math
import sys
from collections import Counter

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# ============================================================
# COUNTRY NORMALIZATION
# ============================================================

REGION_TO_COUNTRY = {
    "SARAWAK": "Malaysia",
    "KCH": "Malaysia",
    "MIRI": "Malaysia",
}

def split_camel_case(s):
    return re.sub(r'(?<!^)(?=[A-Z])', ' ', s)

def normalize_country(raw):
    if not raw or raw.lower() in {"na", ""}:
        return "Unknown"

    raw = raw.split("_")[0].split("/")[0].strip()
    raw = split_camel_case(raw)
    raw_upper = raw.upper()

    if raw_upper in REGION_TO_COUNTRY:
        return REGION_TO_COUNTRY[raw_upper]

    try:
        return pycountry.countries.lookup(raw).name
    except Exception:
        pass

    names = [c.name for c in pycountry.countries]
    best, score = process.extractOne(raw, names)
    return best if score >= 80 else "Unknown"

def country_to_continent(country):
    if country == "Unknown":
        return "Unknown"
    try:
        alpha2 = pycountry.countries.lookup(country).alpha_2
        code = pc.country_alpha2_to_continent_code(alpha2)
        return pc.convert_continent_code_to_continent_name(code)
    except Exception:
        return "Unknown"

# ============================================================
# HEADER PARSER
# ============================================================

YEAR_RE = re.compile(r"(19|20)\d{2}")

def parse_header(header, verbose=False):
    header = header.replace(">", "").strip()
    raw_country = None
    year = None

    if "|" in header:
        for part in header.split("|"):
            if YEAR_RE.fullmatch(part):
                year = int(part)
            elif part.isalpha() or "_" in part:
                raw_country = part

    elif "/" in header:
        raw_country = "Malaysia"
        for part in header.split("/"):
            if part.isdigit():
                year = 2000 + int(part) if len(part) == 2 else int(part)

    country = normalize_country(raw_country)

    if verbose:
        logging.debug(f"{header} → raw={raw_country}, country={country}, year={year}")

    return raw_country or "Unknown", country, year or "Unknown"

def extract_subgenotype(header):
    return header.split("#")[0] if "#" in header else "Unknown"

# ============================================================
# SEQUENCE PROCESSING
# ============================================================

def process_sequence(seq, start=None, stop=None):
    seq = seq.upper().lstrip("N").rstrip("N")

    if start or stop:
        seq = seq[(start - 1 if start else 0):(stop if stop else len(seq))]

    L = len(seq) - (len(seq) % 3)
    seq = seq[:L]

    return seq if len(seq) >= 90 else None

# ============================================================
# KaKs_Calculator RUNNER
# ============================================================

def write_fasta(seq1, seq2, path):
    with open(path, "w") as f:
        f.write(f">seq1\n{seq1}N\n>seq2\n{seq2}N\n")

def run_kaks_worker(args):
    seq1, seq2, gA, gB, genetic_code, method, verbose = args

    with tempfile.NamedTemporaryFile("w+", suffix=".fasta") as fasta_tmp, \
         tempfile.NamedTemporaryFile("w+", suffix=".txt") as out_tmp:

        write_fasta(seq1, seq2, fasta_tmp.name)

        cmd = [
            "KaKs_Calculator",
            "-i", fasta_tmp.name,
            "-o", out_tmp.name,
            "-c", str(genetic_code),
            "-m", method
        ]

        if verbose:
            logging.debug(" ".join(cmd))

        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            if verbose:
                logging.warning(proc.stderr)
            return None

        try:
            df = pd.read_csv(out_tmp.name, sep="\t")
            ka = float(df.loc[0, "Ka"])
            ks = float(df.loc[0, "Ks"])
            return gA, gB, ka, ks
        except Exception:
            return None

# ============================================================
# MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", default="dnds_output.csv")
    ap.add_argument("--gene")
    ap.add_argument("--start", type=int)
    ap.add_argument("--stop", type=int)
    ap.add_argument("--genetic_code", type=int, default=1)
    ap.add_argument("--method", default="YN")
    ap.add_argument("--verbose", action="store_true")

    # ✅ FIX: accept both underscore and hyphen forms
    ap.add_argument(
        "--save_trimmed_fasta",
        "--save-trimmed-fasta",
        dest="save_trimmed_fasta",
        action="store_true",
        help="Save trimmed FASTA sequences used for Ka/Ks analysis"
    )

    args = ap.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    records = []
    country_counts = Counter()

    for rec in SeqIO.parse(args.input, "fasta"):
        raw, country, year = parse_header(rec.description, args.verbose)
        seq = process_sequence(str(rec.seq), args.start, args.stop)
        if not seq:
            continue

        records.append({
            "id": rec.id,
            "seq": seq,
            "country": country,
            "continent": country_to_continent(country),
            "year": year,
            "subgenotype": extract_subgenotype(rec.description)
        })
        country_counts[country] += 1

    logging.info("Countries detected:")
    for c, n in country_counts.items():
        logging.info(f"  {c}: {n}")

    df = pd.DataFrame(records)

    if df.empty:
        logging.error("No usable sequences after filtering")
        sys.exit(1)

    countries = [c for c in df.country.unique() if c != "Unknown"]

    tasks = []

    for c1, c2 in itertools.combinations(countries, 2):
        s1 = df[df.country == c1].seq
        s2 = df[df.country == c2].seq
        for a in s1:
            for b in s2:
                tasks.append((a, b, c1, c2, args.genetic_code, args.method, args.verbose))

    for c in countries:
        seqs = df[df.country == c].seq.tolist()
        for a, b in itertools.combinations(seqs, 2):
            tasks.append((a, b, c, c, args.genetic_code, args.method, args.verbose))

    logging.info(f"Total KaKs tasks: {len(tasks)}")

    with Pool(max(1, cpu_count() - 1)) as pool:
        results = list(tqdm(pool.imap(run_kaks_worker, tasks), total=len(tasks)))

    agg = {}
    for r in results:
        if r:
            agg.setdefault((r[0], r[1]), []).append((r[2], r[3]))

    if not agg:
        logging.error("KaKs_Calculator produced zero usable results")
        sys.exit(1)

    rows = []

    for (gA, gB), vals in agg.items():
        mean_ka = sum(v[0] for v in vals) / len(vals)
        mean_ks = sum(v[1] for v in vals) / len(vals)
        ratio = math.inf if mean_ks == 0 else mean_ka / mean_ks

        base = {
            "Group_A": gA,
            "Group_B": gB,
            "Continent_A": country_to_continent(gA),
            "Continent_B": country_to_continent(gB),
            "mean_Ka": mean_ka,
            "mean_Ks": mean_ks,
            "Ka/Ks": ratio,
            "pairs": len(vals),
            "gene": args.gene,
            "start": args.start,
            "stop": args.stop
        }

        rows.append(base)

        if gA != gB:
            rows.append({
                **base,
                "Group_A": gB,
                "Group_B": gA,
                "Continent_A": base["Continent_B"],
                "Continent_B": base["Continent_A"]
            })

    pd.DataFrame(rows).to_csv(args.output, index=False)
    logging.info(f"Summary output written: {args.output}")

    if args.save_trimmed_fasta:
        out_fasta = os.path.splitext(args.output)[0] + "_trimmed.fasta"
        with open(out_fasta, "w") as f:
            for r in records:
                suffix = r["country"].replace(" ", "_")
                f.write(f">{r['id']}_{suffix}\n{r['seq']}\n")
        logging.info(f"Trimmed FASTA written: {out_fasta}")

if __name__ == "__main__":
    main()

