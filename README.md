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
- `annotation_file`: path to gff3 file with annotation to reference
- `output_prefix`: output folder name and prefix for all pipeline outputs (alignment, annotation and NEXUS file).
