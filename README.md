# MrBayes Preparation Pipeline

This repository contains a simple pipeline that aligns a set of sequences to a reference genome and prepares a partitioned NEXUS file for MrBayes.

## Usage

1. Provide a reference FASTA file and a multi‑FASTA file of sequences to analyse. A GFF3 annotation file with the same base name as the reference is expected alongside it.
2. Edit `config.yaml` to point to your input files and choose an output prefix.
3. Run the pipeline:

```bash
python3 mrbayes_pipeline.py -c config.yaml
```

The script will create an alignment with MAFFT, adjust the reference annotations to match the alignment, fill unannotated regions and finally write a NEXUS file containing a MrBayes block with partitions.

## Configuration file

`config.yaml` contains the following keys:

- `reference_fasta`: path to the reference genome sequence.
- `sequences_fasta`: path to the sequences to align.
- `output_prefix`: prefix for all pipeline outputs (alignment, annotation and NEXUS file).

The reference annotation file is expected at the same location as the reference
FASTA, with the extension replaced by `.gff3`.

Ensure that MAFFT is installed and available on your system path before running the pipeline.

The pipeline preserves `?` characters that represent missing data in the input
FASTA files. It masks them as `N` during the MAFFT run and restores the `?`
positions afterwards so they appear in the final alignment.
