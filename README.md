# resistogrid

`resistogrid` pivots long-format TSV outputs from NCBI AMRFinder and AMRFinderPlus into clean, sample-by-gene matrices.

Many Laboratory Information Management Systems (LIMS) require structured binary matrix inputs (0/1 presence or absence) to process and store Antimicrobial Resistance (AMR) detection results automatically. `resistogrid` bridges the gap between raw bioinformatic tool outputs and LIMS ingestion pipelines by normalizing headers, resolving gene families, and exporting standardized matrices.

---

## Features

* **LIMS-Ready Binary Matrices:** Export matrix cells as binary integer values (`1`/`0`), detection strings (`Detected`/`Not Detected`), `% identity`, `% coverage`, or hit counts.
* **Automatic Schema Detection:** Seamlessly ingests both current AMRFinderPlus outputs and legacy AMRFinder schemas.
* **Gene Family Resolution:** Groups individual alleles into parent resistance families (e.g., mapping `blaOXA-181` to the `blaOXA-48` family) using an embedded, optimized subset of the NCBI Reference Gene Hierarchy.
* **Hierarchical Matrix Formatting:** Exports side-by-side family presence (`blaOXA-48_fam`), list of detected family members (`blaOXA-48_fam_elements`), and individual gene columns.
* **Quality & Subtype Filtering:** Filter hits by minimum coverage, identity, alignment length, element type (e.g., `AMR`, `VIRULENCE`), or detection method (`BLASTX`, `ALLELEX`, etc.).
* **Visualization Outputs:** Generate long-format 3-column TSVs (`Source`, `Family`, `Target`) ready for Sankey diagram generation.

---

## Installation

```bash
pip install resistogrid
```

```bash
conda install -c bioconda resistogrid
```

For local development or installing directly from the repository:

```bash
git clone [https://github.com/your-org/resistogrid.git](https://github.com/your-org/resistogrid.git)
cd resistogrid
pip install -e .
```

## Quick Start

### Basic LIMS Binary Export
Pivot one or more AMRFinder TSV files into a binary matrix (`1` for present, `0` for absent):

```bash
resistogrid -f sample1.tsv sample2.tsv -o lims_amr_matrix.tsv --matrix-value binary
```

## Process an Entire Directory with Quality Thresholds

Recursively scan a folder of TSV files, filter for AMR elements with >= 90% coverage and >= 95% identity, and output a hierarchical matrix:

```bash
resistogrid -f data/amrfinder_results/ \
  -o processed_matrix.tsv \
  --matrix-value binary \
  -t AMR \
  --min-cov 90.0 \
  --min-ident 95.0
```

## Generate Sankey Flow Table

Output a 3-column long table (Source, Family, Target) alongside your pivoted matrix for Sankey diagram rendering:

```bash
resistogrid -f data/*.tsv -o matrix.tsv --sankey-out sankey_links.tsv
```

## Update Embedded NCBI Gene Hierarchy

Update the package's local lookup database directly from the NCBI FTP server:

```bash
resistogrid --update-db
```

## Hierarchical Matrix Output Structure

When using `--group-by hierarchical` (default), `resistogrid` generates a structured matrix format ideal for both automated parsing and human review:

| Sample | blaOXA-48_fam | blaOXA-48_fam_elements | blaOXA-181 | blaCTX-M_fam | blaCTX-M_fam_elements | blaCTX-M-15 |
| :--- | :---: | :--- | :---: | :---: | :--- | :---: |
| Sample_01 | **1** | blaOXA-181 | **1** | **0** | | **0** |
| Sample_02 | **1** | blaOXA-181, blaOXA-232 | **1** | **1** | blaCTX-M-15 | **1** |
| Sample_03 | **0** | | **0** | **1** | blaCTX-M-15 | **1** |

## CLI Reference

```bash
usage: resistogrid [-h] [-v] [--update-db [DIR]] [-f FILES [FILES ...]]
                   [-o OUTPUT] [--sankey-out SANKEY_OUT] [--key-file KEY_FILE]
                   [--matrix-value {detected,binary,identity,coverage,count}]
                   [--group-by {hierarchical,symbol,family,family_symbol}]
                   [-g GENE_FILTER] [-t ELEMENT_TYPE [ELEMENT_TYPE ...]]
                   [-m METHOD [METHOD ...]] [--min-target-len MIN_TARGET_LEN]
                   [--min-ref-len MIN_REF_LEN] [--min-cov MIN_COV]
                   [--min-ident MIN_IDENT] [--min-aln-len MIN_ALN_LEN]
```

### Command Options

**Database Maintenance**
* `--update-db [DIR]`: Download and optimize the latest NCBI gene hierarchy. Defaults to package internal data if `DIR` is omitted.

**Input / Output**
* `-f, --files`: Space-separated list of input TSV files or directory paths.
* `-o, --output`: Output file path for the pivoted matrix (default: `resistogrid_matrix.tsv`).
* `--sankey-out`: Path to write a 3-column TSV (`Source`, `Family`, `Target`) for visualization.
* `--key-file`: Path to a custom gene hierarchy CSV file.

**Matrix Values & Grouping**
* `--matrix-value`: Type of values populating the matrix cells:
  * `detected`: Strings (`"Detected"` / `"Not Detected"`)
  * `binary`: Integers (`1` / `0`)
  * `identity`: Maximum % Identity float value (`0.0` if absent)
  * `coverage`: Maximum % Coverage float value (`0.0` if absent)
  * `count`: Total hit count integer
* `--group-by`: Matrix column organization:
  * `hierarchical`: Creates family summary, element list, and gene columns
  * `symbol`: Direct gene symbol columns (`blaOXA-181`)
  * `family`: Aggregated gene family columns (`blaOXA-48`)
  * `family_symbol`: Prefixed columns (`blaOXA-48___blaOXA-181`)

**Filtering**
* `-g, --gene-filter`: Case-insensitive substring search matching Gene Symbol or Sequence Name.
* `-t, --element-type`: Filter by element type (`AMR`, `VIRULENCE`, `STRESS`, `POINT`).
* `-m, --method`: Filter by identification method (`BLASTX`, `EXACTX`, `ALLELEX`, etc.).
* `--min-cov`: Minimum % Coverage threshold (`0.0` - `100.0`).
* `--min-ident`: Minimum % Identity threshold (`0.0` - `100.0`).
* `--min-target-len`: Minimum target sequence length.
* `--min-ref-len`: Minimum reference sequence length.
* `--min-aln-len`: Minimum alignment length.

## Python API Usage

You can also call resistogrid functions directly within your Python scripts:

```python
from pathlib import Path
from resistogrid.converter import process_amrfinder_files

input_files = [Path("data/sample1.tsv"), Path("data/sample2.tsv")]
output_matrix = Path("output/lims_matrix.tsv")

process_amrfinder_files(
    files=input_files,
    output=output_matrix,
    matrix_value="binary",
    group_by="hierarchical",
    min_cov=90.0,
    min_ident=95.0,
)
```

## AI Assistance & Verification Disclaimer

This software was developed with the assistance of artificial intelligence (AI) coding tools. While the codebase includes automated unit tests and validation procedures, bioinformatic outputs processed by `resistogrid` are provided on an **"AS IS" basis, without warranties or conditions of any kind**.

Users—especially those in clinical, diagnostic, or public health laboratory environments—are responsible for independently validating all matrix outputs, gene family mappings, and LIMS integrations prior to making clinical, diagnostic, or regulatory decisions.

## License
Distributed under the GNU General Public License v3 (GPLv3).