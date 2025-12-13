# Pairwise dN/dS Analysis Script

## About

This Python program calculates pairwise **dN (nonsynonymous substitutions), dS (synonymous substitutions), and dN/dS ratios** between nucleotide coding sequences in a FASTA file. It allows for grouping sequences by country or country-year, trimming sequences, and optional gene start/stop positions.  

Key features:

* Supports FASTA files with variable header formats
* Normalizes country names automatically using `pycountry`
* Optionally trims sequences of leading/trailing N's
* Supports pairwise comparisons grouped by **country** or **country-year**
* Uses Biopython codon alignment (`Bio.codonalign`) to calculate dN, dS, and dN/dS without relying on external `codeml`
* Multiprocessing enabled for efficient computation
* Verbose mode for debugging and diagnostic output

---

## Installation Instructions

### Using Conda (Recommended)

1. Ensure you have **Conda** installed (Anaconda or Miniconda)
2. Create a new environment using the provided YAML file:
```bash
conda env create -f environment.yml
conda activate dn_ds_analysis
```

## This installs:
- Python 3.11
- Biopython
- pandas
- pycountry
- thefuzz
- tqdm
Optional: python-Levenshtein for faster fuzzy matching with thefuzz.
```bash
### Required Arguments:
- -i, --input: Path to input FASTA file containing coding sequences
- -o, --output: Output CSV file for pairwise dN, dS, and dN/dS values

### Optional Arguments:
- --gene: Gene name for metadata in output CSV
- --start: Start position (1-based) to trim sequences before analysis
- --stop: Stop position (1-based) to trim sequences before analysis
- --mode: Grouping mode for comparisons:
    - country (default) — compare sequences by country
    - country-year — compare sequences by country-year combination
- --verbose: Enable verbose/debug logging
- --threads: Number of CPU threads to use (default: all available minus one)
```
## Output
The script produces:
1. Pairwise CSV (OUTPUT_CSV) with columns:
*Group_A*	*Group_B*	*mean_dN*	*mean_dS*	*dN/dS*	*pairs*	*gene*	*start*	*stop*
  - Group_A and Group_B correspond to country or country-year depending on mode
  - mean_dN / mean_dS are averages over all site comparisons between sequence pairs
  - dN/dS is the ratio of mean dN to mean dS
  - pairs is the number of pairwise comparisons aggregated
2. Verbose logging prints detected countries, sequence counts, and diagnostic info


## Example
Given a FASTA file sequences.fasta:
```bash
>D5/8242|CVA6|SWK/OCT/08
ATGGGCGCCCAAGTTTCAACAGAAAAGTCTGGGTCGCACGAGACAAAGAATGTAGCGACTGAAGGGTCTACTATCAACTTCACCAACATCAATTACTATAA
>D3/KR815992|CVA6|China|2013
ATGGGCGCCCAAGTCTCAACAGAAAAATCTGGGTCGCACGAGACAAAGAATGTAGCGACTGAAGGGTCTACTATCAACTTCACCAACATCAATTACTATAA
>D3/MZ546175|CVA6|China|2016
ATGGGCGCCCAAGTTTCAACAGAAAAATCTGGGTCGCACGAGACAAAGAATGTAGCGACCGAAGGGTCTACTATCAACTTCACCAACATCAATTACTATAA
```
Run the analysis:
```bash
python dn_ds_analysis.py -i sequences.fasta -o dnds_output.csv --mode country --verbose --gene VP1
```

## Notes
--Ensure all input sequences are coding sequences with length multiple of 3 for codon alignment
--Leading/trailing Ns are trimmed automatically
--Multiprocessing will use all available cores minus one by default
--Verbose mode is useful to check country normalization and intermediate pairwise results
--Gene, start, and stop positions can be used to trim the sequence before analysis for consistent region comparisons
