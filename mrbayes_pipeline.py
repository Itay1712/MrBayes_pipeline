import yaml
from pathlib import Path
from Bio import AlignIO
from annotate_sequences import (
    run_mafft_alignment_linux,
    remove_gaps_and_ns_from_alignment,
    adjust_annotation_for_gaps_gff3,
    fill_unannotated_regions_gff3,
)
from nexus_from_annotation import create_nexus_from_alignment_and_annotation


def run_pipeline(config_file: str):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    ref_fasta = cfg['reference_fasta']
    seqs_fasta = cfg['sequences_fasta']
    annotation = str(Path(ref_fasta).with_suffix('.gff3'))
    prefix = cfg.get('output_prefix', 'output')
    include_partitions = cfg.get('include_partitions')

    alignment_file = f"{prefix}_alignment.fasta"
    clean_alignment = f"{prefix}_clean.fasta"
    adjusted_annotation = f"{prefix}_adjusted.gff3"
    filled_annotation = f"{prefix}_filled.gff3"
    nexus_file = f"{prefix}.nex"

    # Align sequences to reference
    run_mafft_alignment_linux(ref_fasta, seqs_fasta, alignment_file)
    remove_gaps_and_ns_from_alignment(alignment_file, clean_alignment)

    # Adjust annotation to include gaps from alignment
    adjust_annotation_for_gaps_gff3(annotation, clean_alignment, adjusted_annotation)

    # Calculate alignment length and fill unannotated regions
    alignment = AlignIO.read(clean_alignment, 'fasta')
    alignment_length = alignment.get_alignment_length()
    fill_unannotated_regions_gff3(adjusted_annotation, alignment_length, filled_annotation)

    # Create NEXUS file with partitions for MrBayes
    create_nexus_from_alignment_and_annotation(
        clean_alignment,
        filled_annotation,
        nexus_file,
        include_partitions=include_partitions,
    )


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run MrBayes preparation pipeline')
    parser.add_argument('-c', '--config', default='config.yaml', help='YAML configuration file')
    args = parser.parse_args()

    run_pipeline(args.config)
