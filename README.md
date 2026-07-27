# Bacterial Outbreak Detection Pipeline

A Nextflow pipeline that takes paired-end bacterial sequencing data (from NCBI SRA and/or local FASTQs), assembles and QC's each genome, identifies each sample's genus/species, genotypes it, and builds a core-genome SNP phylogenetic tree per genus to relate samples to each other for outbreak detection.

```
                                    ┌─→ fastqMetrics
sra.txt ──→ downloadSRA ──→ convertSRA ─┐
                                          ├─→ cleanFastq ──→ assembleGenome ──→ filterContigs ──┬─→ genePrediction
input/*_R{1,2}_*.fastq.gz ───────────────┘                                                       │
                                                                                                   ├─→ [optional 16S branch, off by default:
                                                                                                   │      extract16S → genusIdentification (BLAST, genus only)
                                                                                                   │        → genusRefList (NCBI datasets) → speciesIdentification (FastANI, species-level)
                                                                                                   │      every sample's 16S sequence, regardless of genus, also feeds a shared tree:
                                                                                                   │        combine16SSequences → align16SSequences (MAFFT) → build16STree (IQ-TREE2) → render16STree]
                                                                                                   │
                                                                                                   ├─→ downloadMashDB → mashIdentification → readMashIdentification (genus + species)
                                                                                                   │        ├─→ genotyping (mlst)
                                                                                                   │        ├─→ genomeQA (BUSCO)
                                                                                                   │        └─→ contigQA (GUNC)
                                                                                                   └─→ grouped by genus → phylogeneticAnalysis (Parsnp: alignment + core genome)
                                                                                                                            → buildCoreGenomeTree (IQ-TREE2) → renderPhylogeneticTree (FigTree/SVG)
                                                                                                                                  │
                                                                                                            all per-sample results ─→ summarizeResults → combineSummaryReports
```

Reads can be supplied as an SRA accession list, a directory of local paired-end FASTQs, or both at once.

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| Nextflow | 25.10.4 | Latest release |
| Miniforge (conda) | 26.1.0 |Package manager |
| Java | 23.0.2 | Required by Nextflow |
| Docker Desktop | latest | Required for the Parsnp phylogenetics step (see below) |
| OS | macOS (Apple Silicon arm64) | Tested on macOS Sequoia 15 |

> Most pipeline tool dependencies (fastp, skesa, seqkit, sra-tools, entrez-direct, pigz, prodigal, barrnap, bedtools, blast+, mash, mlst, busco, gunc, mafft, iqtree2, figtree, ncbi-datasets-cli, fastani, and their Python dependencies) are installed automatically via conda using `nf_cmds.config` — no manual tool installation needed beyond the requirements above.
>
> **Exception: Parsnp.** Parsnp has no native osx-arm64 conda build, so the `phylogeneticAnalysis` process runs it inside a Docker container (`staphb/parsnp:1.5.6`, linux/amd64, transparently emulated by Docker Desktop on Apple Silicon) instead of via conda. Docker Desktop must be installed and running for this step to work. Everything downstream including tree rendering (`renderPhylogeneticTree`, via FigTree) runs natively through conda.

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

**4. Install Docker Desktop (required for the Parsnp step)**

Download and install from [docker.com](https://www.docker.com/products/docker-desktop/), then make sure it's running before launching the pipeline. The `staphb/parsnp:1.5.6` image will be pulled automatically on first use.

---

## Workflow Test

This pipeline accepts reads from two independent sources — an SRA accession list and/or a local directory of paired-end FASTQs — which merge into the same run. The test dataset below exercises both pathways at once.

### Test Dataset

**SRA pathway:** a text file input listing SRA IDs. The test data, **sra.txt**, includes two small bacterial accessions: **SRR1172848** (*Mycobacterium tuberculosis H3367*) and **SRR2584863** (*Escherichia coli B str. REL606*).

**Local-directory pathway:** a directory of local paired-end FASTQs, **input/**, containing two samples' reads, matching the naming convention: `<sample>_S##_L###_R{1,2}_001.fastq.gz`, e.g.:
```
input/
├── Aae52ace0b_S01_L001_R1_001.fastq.gz
├── Aae52ace0b_S01_L001_R2_001.fastq.gz
├── Af62ef3a41_L001_R1_001.fastq.gz
└── Af62ef3a41_L001_R2_001.fastq.gz
```


Note that across both pathways this test set spans multiple genera. Since `phylogeneticAnalysis` groups samples by their Mash-identified genus before alignment, and both Parsnp's alignment and the downstream IQ-TREE2 tree need at least 3 genomes for a real topology, any genus with fewer than 3 samples in your test set will produce a placeholder "no tree" result rather than an actual phylogeny. Use the pipeline on 3+ same-genus samples to see a real core-genome tree.

### Test Commands
Clone or download this repository and navigate into it, then run the following **2 commands**:

**Step 1 — Run the pipeline**
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt
```

**Step 2 — Confirm outputs**
```bash
tree -L 2 results/raw_fastq results/clean_fastq results/Assemblies results/summary_reports
```

If you only want to test one pathway at a time, either flag can be run alone:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt        # SRA only
nextflow run nf_cmds.nf -c nf_cmds.config --input_dir ./input        # local FASTQs only
```

---

## Running on Your Own Data

**From SRA accessions:** create a plain text file with one accession per line:
```bash
printf "SRR2584863\nSRR9094324\nSRR1172848\nSRR2093876\n" > sra.txt
```
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt
```

**From local FASTQs:** point `--input_dir` at a directory of paired-end reads named like `Sample_S01_L001_R1_001.fastq.gz` / `..._R2_001.fastq.gz`:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --input_dir ./input
```

**Both at once** — SRA downloads and local FASTQs are merged into the same run:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt --input_dir ./input
```

For meaningful outbreak/relatedness trees, include at least 3 samples per genus of interest; Parsnp's alignment step and the downstream IQ-TREE2 tree need a minimum of 3 genomes to produce a real topology .

### Optional: 16S rRNA identification

By default, genus/species identification is done via Mash which is a whole-assembly comparison against a local NCBI RefSeq sketch database. This prevents failures from fragmented assemblies typical of short-read sequencing. A 16S rRNA-based identification side-branch is available for but is off by default, since 16S is fragile to recover from fragmented short-read assemblies and multi-copy rRNA operons often don't assemble cleanly. It does **not** feed into genotyping or the final summary report either way. Turning it on runs three additional stages:

1. **Genus ID** — `extract16S` (barrnap) pulls each sample's 16S sequence, then `genusIdentification` BLASTs it against NCBI's 16S rRNA database. BLASTn alone only reliably resolves down to genus.
2. **Species confirmation** — `genusRefList` fetches reference genomes for that genus via the NCBI datasets CLI, and `speciesIdentification` runs whole-genome FastANI against them to resolve species (finer resolution than 16S BLAST alone can offer).
3. **A shared 16S tree** — every sample's 16S sequence, regardless of genus, is combined (`combine16SSequences`), aligned with MAFFT (`align16SSequences`), and built into a single IQ-TREE2 tree (`build16STree` → `render16STree`). Unlike the per-genus core-genome tree, this one is meant to span diverse taxa; it's a broad relatedness check across your whole sample set, not an outbreak-resolution tree (ie a complement, not a replacement).

Every step in this branch has `errorStrategy = 'ignore'` set, so a failure anywhere in the chain (ie no 16S hit for one sample) only drops that branch for that sample rather than failing the whole run. Steps needing live NCBI access (`genusIdentification`, `genusRefList`, `readSpeciesIdentification`) require internet connectivity.

To turn it on:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config --sra_list sra.txt --run_16S_id true
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--sra_list` | `./sra.txt` | Text file of SRA accessions, one per line |
| `--input_dir` | `./input` | Directory of local paired-end `*_R{1,2}_*.fastq.gz` reads |
| `--mash_db_dir` | `${projectDir}/data/mash` | Where the ~950 MB Mash RefSeq sketch + assembly summary are cached (downloaded once, reused across runs) |
| `--run_16S_id` | `false` | Turn on the optional BLAST-based 16S rRNA genus ID side-branch |

At least one of `--sra_list` or `--input_dir` must point to real input, or the pipeline will exit with an error.

---

## Output Files

| Directory | File(s) | Description |
|-----------|---------|-------------|
| `raw_fastq/` | `*_1.fastq.gz`, `*_2.fastq.gz` | Raw paired-end reads from SRA |
| `clean_fastq/` | `*.R1.fq.gz`, `*.R2.fq.gz` | Quality-trimmed reads |
| `clean_fastq/` | `*.json`, `*.html` | fastp QC reports |
| `clean_fastq/` | `*_stats.tsv` | seqkit read statistics (N50, Q30, etc.) |
| `Assemblies/` | `*.fna` | Raw SKESA genome assemblies |
| `filtered_assemblies/` | `*.fna`, `*_discarded_contigs.fna`, `*.log` | Assemblies with short (<1000 bp) and low-coverage (<10X) contigs removed |
| `gene_prediction/` | `*.gff`, `*.log` | Prodigal ab initio gene predictions |
| `16S_rRNA/` | `*_16S.gff`, `*_16S.fa`, `*_16S_rRNA_blastn.tsv` | Optional: barrnap 16S coordinates/sequence and BLAST genus identification (only if `--run_16S_id true`) |
| `genus_ref_list/` | `*_ref_genomes/*.fna` | Optional: reference genomes fetched per genus for FastANI (only if `--run_16S_id true`) |
| `species_identification/` | `*-fastani.tsv` | Optional: FastANI species confirmation against genus reference genomes (only if `--run_16S_id true`) |
| `16S_tree/` | `16S_tree.treefile`, `16S_tree.svg` | Optional: single IQ-TREE2 tree across every sample's 16S sequence, spanning all genera (only if `--run_16S_id true`) |
| `mash_identification/` | `*_mash_taxonomy_results.tsv` | Mash-based genus/species identification against NCBI RefSeq |
| `genotyping_results/` | `*-genotyping.tsv` | MLST scheme/sequence type per sample |
| `genomeQA/` | `*-busco_summary.txt` | BUSCO genome completeness assessment |
| `contigQA/` | `*.tsv` | GUNC contamination/chimerism assessment |
| `phylogenetic_analysis/` | `*-parsnp.tree`, `*-core_alignment.fasta`, `*-core_genome.treefile`, `*-core_genome.svg` | Per-genus core-genome workflow: Parsnp produces the alignment (and its own quick tree, kept for reference), IQ-TREE2 builds the tree actually used downstream, FigTree renders it to SVG |
| `summary_reports/` | `summary_report.tsv` | One merged report row per sample: genotyping, QC, QA, taxonomy, and phylogenetic tree |

---

## Resuming a Failed Run

Nextflow caches completed tasks. If the pipeline fails mid-run, resume from where it stopped with:
```bash
nextflow run nf_cmds.nf -c nf_cmds.config -resume
```

## Citation
Created for GA Tech BIOL7210
