# Bacterial Genome Assembly Pipeline

A Nextflow pipeline that downloads paired-end bacterial sequencing data from NCBI SRA, performs quality trimming with fastp, assembles genomes with SKESA, and collects read metrics with seqkit. This pipeline illustrates the **sequential and parallel execution** of Nextflow.

```
SRA Accession List → downloadSRA → convertSRA → cleanFastq → assembleGenome
                                                     ↓
                                               fastqMetrics
```

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| Nextflow | 25.10.4 | Latest release |
| Miniforge (conda) | 26.1.0 | Recommended package manager |
| Java | 23.0.2 | Required by Nextflow |
| OS | macOS (Apple Silicon arm64) | Tested on macOS Sequoia 15 |

> All pipeline tool dependencies (fastp, skesa, seqkit, sra-tools, pigz) are installed automatically via conda using `nf_cmds.config` (no manual tool installation needed beyond the requirements above).

---

## Installation

**1. Install Java (if not already installed)**
```bash
sdk install java 21-tem
```
> If you don't have SDKMAN: `curl -s "https://get.sdkman.io" | bash`

**2. Install Nextflow**
```bash
curl -s https://get.nextflow.io | bash
mkdir -p $HOME/.local/bin && mv nextflow $HOME/.local/bin/
echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.zshrc && source ~/.zshrc
```

**3. Install Miniforge (if not already installed)**
```bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh -o miniforge.sh
bash miniforge.sh
```

---

## Workflow Test

### Test Dataset
This pipeline takes a text file input that includes a list of SRA IDs. The test data, **sra.tx**t, includes two small bacterial accessions: **SRR1172848** (*Mycobacterium tuberculosis H3367*) and **SRR2584863** (*Escherichia coli B str. REL606*). This test data should take less than 20 minutes to run and less than 0.4 CPU hours. 

### Test Commands 
Clone or download this repository and navigate into it, then run the following **2 commands**:

**Step 1 — Run the pipeline**
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt
```

**Step 2 — Confirm outputs**
```bash
tree -L 2 results/raw_fastq results/clean_fastq results/Assemblies
```

Expected output structure:
```
results/raw_fastq/
├── SRR2584863_1.fastq.gz
├── SRR2584863_2.fastq.gz
├── SRR1172848_1.fastq.gz
└── SRR1172848_2.fastq.gz
results/clean_fastq/
├── SRR2584863.R1.fq.gz
├── SRR2584863.R2.fq.gz
├── SRR2584863.json
├── SRR2584863.html
├── SRR2584863_stats.tsv
├── SRR1172848.R1.fq.gz
├── SRR1172848.R2.fq.gz
├── SRR1172848.json
├── SRR1172848.html
└── SRR1172848_stats.tsv
results/Assemblies/
├── SRR2584863.fna
└── SRR1172848.fna
```

---

## Running on Your Own Data

Create a plain text file with one SRA accession per line:
```bash
printf "SRR2584863\nSRR9094324\nSRR1172848\nSRR2093876\n" > sra.txt
```
Example:
```
SRR2584863
SRR9094324
SRR1172848
SRR2093876
```

Then pass it to the pipeline with `--sra_list`:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt
```

---

## Output Files

| Directory | File(s) | Description |
|-----------|---------|-------------|
| `raw_fastq/` | `*_1.fastq.gz`, `*_2.fastq.gz` | Raw paired-end reads from SRA |
| `clean_fastq/` | `*.R1.fq.gz`, `*.R2.fq.gz` | Quality-trimmed reads |
| `clean_fastq/` | `*.json`, `*.html` | fastp QC reports |
| `clean_fastq/` | `*_stats.tsv` | seqkit read statistics (N50, Q30, etc.) |
| `Assemblies/` | `*.fna` | Assembled genome contigs (≥ 1000 bp) |

---

## Resuming a Failed Run

Nextflow caches completed tasks. If the pipeline fails mid-run, resume from where it stopped with:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config -resume
```

## Citation
Created for GA Tech BIOL7210
