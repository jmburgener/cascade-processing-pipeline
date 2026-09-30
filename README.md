# Cascade Processing Pipeline

[![CI](https://github.com/jmburgener/cascade-processing-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/jmburgener/cascade-processing-pipeline/actions/workflows/ci.yml)

A Nextflow pipeline that takes paired-end, UMI-tagged bisulfite sequencing libraries from FASTQ
to deduplicated alignments, per-CpG methylation calls, and a cohort-level QC summary, with a
spike-in-based QC gate that drops failed libraries before summarizing.

Coding sample on synthetic data, not a production tool. Everything, including the reference,
runs from the repo with no external data.

## Workflow

```
FASTQ (R1/R2)
  │
  ├─ extractUMI          umi_tools extract: move the inline UMI into the read name
  ├─ trimAdapters        cutadapt: 3' adapter (+ UMI), quality trim, length filter
  ├─ bismarkAlign        Bismark: align to the bisulfite-converted synthetic reference
  ├─ sortBam             samtools sort + index
  ├─ dedupUMI            umi_tools dedup: collapse PCR duplicates by position + UMI
  ├─ nameSortBam         samtools sort -n: mates adjacent for methylation calling
  ├─ methylationExtract  bismark_methylation_extractor: per-CpG coverage file
  ├─ conversionQC        conversion efficiency on an unmethylated spike-in
  │     └─ QC gate       samples below min_conversion (95%) are dropped, with a warning
  └─ sampleSummary       on-target rate + CpG methylation over the panel
        └─ cohort_summary.tsv
```

`bismarkGenomePrep` builds the bisulfite-converted index once and shares it across samples.

**Why a spike-in gate:** the spike-in region is unmethylated by design, so any methylation called
there is unconverted cytosine. `100 − % methylation` in the spike-in estimates bisulfite
conversion efficiency per sample. A library that converted poorly would otherwise look
*hyper*methylated, which is a false positive for a methylation-based assay.

Every process runs in a pinned BioContainers image (Bismark 3.1.0, umi_tools 1.1.6,
cutadapt 5.2, samtools 1.24).

## Requirements

- Nextflow (tested with 26.04.6; uses typed `params` and workflow `output` blocks)
- Docker
- Python 3 (standard library only), only to regenerate the synthetic data

## Quick start

```bash
git clone https://github.com/jmburgener/cascade-processing-pipeline
cd cascade-processing-pipeline
nextflow run scripts/main.nf
```

Run from the repo root: parameter paths are relative to it. Outputs go to `results/`.

## Output

- `results/bams/<sample>.dedup.bam`: deduplicated alignments
- `results/cohort_summary.tsv`: one row per sample that passed QC

With the bundled synthetic cohort, `qc_fail` (simulated at 40% conversion) is dropped at the gate:

```
WARN: qc_fail failed QC, conversion efficiency 45.00% < 95.0%
```

```
sample         conversion_pct  dedup_pairs  on_target_pairs  on_target_frac  on_target_cpgs  on_target_meth_pct
high_signal_1  98.33           14           12               0.857           103             81.19
high_signal_2  95.00           14           12               0.857           104             74.48
low_signal_1   98.33           8            6                0.750           64              21.62
low_signal_2   100.00          8            6                0.750           96              21.38
low_signal_3   100.00          8            6                0.750           77              17.22
```

`high_signal_2` sits exactly at the 95% threshold and passes: the gate is inclusive (`>=`).

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `samplesheet` | `data/nextflow_sample_sheet.txt` | CSV: `sample,fastq_1,fastq_2` |
| `fasta_dir` | `data` | Directory holding the reference FASTA |
| `adapter` | `data/adapter.txt` | FASTA with the 3' adapter sequence |
| `panel` | `data/synthetic_panel.bed` | Target regions for on-target and methylation summaries |
| `spikein` | `data/synthetic_spikein.bed` | Unmethylated spike-in region used for conversion QC |
| `umi_pattern` | `NNX` | umi_tools `--bc-pattern` for both reads |
| `min_conversion` | `95.0` | Minimum conversion efficiency (%) to pass QC |

Override any of them on the command line, e.g. `--min_conversion 98`.

## Synthetic data

`scripts/simulate_data.py` generates the synthetic reference, the target panel, the adapter
file and a small cohort of paired FASTQs (high- and low-methylation samples, plus one
poorly converted library), with UMIs and adapters added and a fixed seed (77) for determinism.
The generated files are committed, so this is only needed to regenerate them:

```bash
python3 scripts/simulate_data.py
```

## Layout

```
scripts/main.nf            pipeline
scripts/simulate_data.py   synthetic data generator
nextflow.config            Docker on
data/                      synthetic reference, panel, spike-in, sample sheets, FASTQs
```

## Known limitations

- Tiny synthetic cohort: the numbers show the pipeline works, not assay performance.
- The deduplicated BAM of a sample that fails QC is still published to `results/bams/`. Only
  the summary excludes it.
- `bismarkGenomePrep` writes the converted index next to the input FASTA (`data/Bisulfite_Genome/`,
  gitignored) rather than to its own work directory.
- CI runs the whole pipeline on every push and checks the QC gate (5 samples pass, `qc_fail`
  is excluded). There are no per-process unit tests (e.g. nf-test) yet.
