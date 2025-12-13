# Pairwise dN/dS (Ka/Ks) Analysis Script
## About

This Python program calculates pairwise **dN (nonsynonymous substitutions), dS (synonymous substitutions), and dN/dS ratios** between nucleotide coding sequences in a FASTA file using KaKs_Calculator

.

## Key features:

* Supports FASTA files with variable header formats, including missing or malformed country/year information.
* Automatically normalizes country names using pycountry and detects continent information.
* Extracts subgenotype from sequence names (prefix before # in headers, e.g., D3a#KP289365 → subgenotype D3a).
* Optionally trims sequences based on start/stop positions and removes leading/trailing Ns.
* Supports between-country and within-country pairwise comparisons.
* Aggregates results as averages for all pairwise sequence comparisons between groups.
* Optional output of trimmed FASTA sequences with renamed headers: existing sequence ID suffixed by simplified country name.
*  Multiprocessing enabled for efficient computation.
*  Verbose mode for detailed logging and diagnostics.

## Installation Instructions

### Using Conda (Recommended)

1. Ensure you have Conda installed (Anaconda or Miniconda).

2. Create a new environment:

```bash
conda env create -f environment.yml
conda activate dn_ds_analysis
```
### This installs:
	- Python 3.11
	- Biopython
	- pandas
	- pycountry
	- pycountry_convert
	- thefuzz (optional: with python-Levenshtein for faster fuzzy matching)
	- tqdm

Make sure KaKs_Calculator is installed and available in your system path.

Usage
```bash
python dn_ds_analysis.py -i INPUT_FASTA -o OUTPUT_CSV [options]
```

### Required Arguments

 - -i, --input — Path to input FASTA file containing coding sequences.

 - -o, --output — Output CSV file for aggregated pairwise dN, dS, and dN/dS values.

### Optional Arguments

 - --gene — Gene name to annotate output CSV.

 - --start — Start position (1-based) to trim sequences before analysis.

 - --stop — Stop position (1-based) to trim sequences before analysis.

 - --verbose — Enable verbose logging for debugging and diagnostics.

 - --threads — Number of CPU threads for multiprocessing (default: all available minus one).

 - --save_trimmed_fasta — Save trimmed sequences in a FASTA file with headers suffixed by simplified country name.

## Input Requirements

* Input sequences must be coding sequences (CDS).

* Sequence lengths should preferably be multiples of 3 (codon-aligned).

* Leading/trailing Ns are automatically trimmed.

Run the analysis:
```bash
python dn_ds_analysis.py -i sequences.fasta -o dnds_output.csv --mode country --verbose --gene VP1
```

Country names are normalized automatically; malformed or missing countries are set to Unknown.

Continent information is inferred automatically from the country.

Subgenotype is extracted from sequence names as the prefix before #.

## Output
### 1. Pairwise CSV (OUTPUT_CSV)


|**Column**	|	**Description**							                        |
|-----------|------------------------------------------------------------------ |
|Group_A	|Name of first comparison group (country)				            |
|Group_B	|Name of second comparison group (country)				            |
|Continent_A|Continent of Group_A					    	                	|
|Continent_B|Continent of Group_B						    	                |
|mean_Ka	|Average dN across all pairwise sequence comparisons between groups	|
|mean_Ks	|Average dS across all pairwise sequence comparisons between groups	|
|Ka/Ks		|Ratio of mean_Ka to mean_Ks						                |
|pairs		|Number of pairwise comparisons aggregated				            |
|gene		|Gene name (optional)							                    |
|start		|Start position used for trimming (optional)stop			        |				
|Stop 		|position used for trimming (optional)				                |

* Between-country comparisons: Every sequence from Group_A compared with every sequence from Group_B.

* Within-country comparisons: Every sequence is compared with every other sequence from the same country.

* Results are aggregated to average Ka/Ks for each group pair.

### 2. Optional Trimmed FASTA

If --save_trimmed_fasta is enabled, a FASTA file is saved with sequences:

Trimmed according to start and stop positions.

Headers renamed with suffix of simplified country name, e.g.:
```bash
>D3a#KR815992_China
ATGGGCGCCCAAGTTTCAACAGAAAAGTCT...
>D3a#MZ546175_China
ATGGGCGCCCAAGTTTCAACAGAAAAGTCT...
```
## Example

Input FASTA sequences.fasta:
```bash
>D5/8242|CVA6|SWK/OCT/08
ATGGGCGCCCAAGTTTCAACAGAAAAGTCTGGGTC...
>D3/KR815992|CVA6|China|2013
ATGGGCGCCCAAGTCTCAACAGAAAAATCTGGGTC...
>D3/MZ546175|CVA6|China|2016
ATGGGCGCCCAAGTTTCAACAGAAAAATCTGGGTC...
```

Run analysis:
```bash
python dn_ds_analysis.py \
    -i sequences.fasta \
    -o dnds_output.csv \
    --gene VP1 \
    --start 1 --stop 612 \
    --mode country \
    --verbose \
    --save_trimmed_fasta
```

## Notes
* Multiprocessing uses all available CPU cores minus one by default.
* Verbose mode prints detected countries, subgenotypes, continents, and sequence counts.
* Trimmed sequences in FASTA include simplified country suffix.
* Ka/Ks is calculated via KaKs_Calculator; sequences are padded with N internally to avoid parsing errors.
* Ensure all input sequences are CDS (length multiple of 3) for accurate Ka/Ks estimation.
