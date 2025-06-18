import subprocess
import argparse
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os


def run_mafft_alignment_linux(reference_fasta, sequences_fasta, output_alignment_file):
    """Run MAFFT while preserving ``?`` characters in the sequences."""

    # Read reference and sequences
    records = list(SeqIO.parse(reference_fasta, "fasta"))
    records += list(SeqIO.parse(sequences_fasta, "fasta"))

    # Record ``?`` positions per sequence and replace them with ``N``
    qmark_positions = {}
    for rec in records:
        seq_str = str(rec.seq)
        positions = [i for i, c in enumerate(seq_str) if c == "?"]
        if positions:
            qmark_positions[rec.id] = set(positions)
            seq_str = seq_str.replace("?", "N")
            rec.seq = Seq(seq_str)

    all_sequences_file = "processed/combined.fasta"
    SeqIO.write(records, all_sequences_file, "fasta")

    mafft_command = f"mafft --auto {all_sequences_file} > {output_alignment_file}"
    subprocess.run(mafft_command, shell=True, check=True)
    print(f"MAFFT alignment completed. Output saved to {output_alignment_file}")
    if os.path.exists(all_sequences_file):
        os.remove(all_sequences_file)

    # Restore ``?`` characters in the aligned sequences if necessary
    if qmark_positions:
        alignment = AlignIO.read(output_alignment_file, "fasta")
        restored_records = []
        for rec in alignment:
            positions = qmark_positions.get(rec.id)
            if not positions:
                restored_records.append(rec)
                continue
            seq_chars = list(str(rec.seq))
            ungapped_idx = 0
            for i, ch in enumerate(seq_chars):
                if ch != "-":
                    if ungapped_idx in positions:
                        seq_chars[i] = "?"
                    ungapped_idx += 1
            restored_records.append(SeqRecord(Seq("".join(seq_chars)), id=rec.id, description=""))

        AlignIO.write(MultipleSeqAlignment(restored_records), output_alignment_file, "fasta")
        print("Restored '?' characters in alignment.")

def remove_gaps_and_ns_from_alignment(input_alignment_file, output_alignment_file):
    """
    Removes columns in the alignment that contain only gaps ('-') and 'N' from a FASTA alignment file.
    
    Parameters:
    - input_alignment_file: Path to the input alignment file in FASTA format.
    - output_alignment_file: Path to save the cleaned alignment in FASTA format.
    """

    alignment = AlignIO.read(input_alignment_file, "fasta")
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    ## collect the columns that are not empty
    clean_columns = []
    for i in range(alignment_length):
        column = [alignment[j, i] for j in range(num_sequences)]

        if not all(c in ['-', 'N', 'n', '?'] for c in column):
            clean_columns.append(i)
    
    cleaned_alignment = []
    for record in alignment:
        cleaned_seq = ''.join([record.seq[i] for i in clean_columns])
        cleaned_alignment.append(SeqRecord(Seq(cleaned_seq), id=record.id, description=""))
    
    AlignIO.write(MultipleSeqAlignment(cleaned_alignment), output_alignment_file, "fasta")
    print(f"Cleaned alignment written to {output_alignment_file}")


def adjust_annotation_for_gaps_gff3(annotation_file, alignment_file, output_file):
    """
    Adjusts GFF3 annotations based on gaps introduced in the alignment.
    
    Parameters:
    - annotation_file: Path to the reference GFF3 annotation file.
    - alignment_file: Path to the multiple sequence alignment (in FASTA format).
    - output_file: Path to the output file with adjusted GFF3 annotations.
    """
    
    alignment = AlignIO.read(alignment_file, "fasta")
    reference_seq = str(alignment[0].seq)
    adjusted_annotations = []

    with open(annotation_file, 'r') as ann_file:
        for line in ann_file:
            if line.startswith("#") or len(line.strip()) == 0:
                adjusted_annotations.append(line.strip())
                continue
            
            fields = line.strip().split("\t")
            
            region_start = int(fields[3])
            region_end = int(fields[4])
                        
            adjusted_start = add_gaps_before(reference_seq, region_start)
            adjusted_end =  add_gaps_before(reference_seq, region_end)

            fields[3] = str(adjusted_start)
            fields[4] = str(adjusted_end)
            
            feature_type = fields[2]

            if feature_type in ["CDS", "exon", "region"]:
                continue

            if feature_type == "tRNA":
                attributes = fields[8]
                product = extract_attribute(attributes, "product")
                codons =  extract_attribute(attributes, "codons")
                if product:
                    fields[2] = f"{product}"
                if codons:
                    fields[2] += f"_{codons}"


            if feature_type == "rRNA":
                attributes = fields[8]
                product = extract_attribute(attributes, "product")
                if product:
                    fields[2] = f"{product}"
            
            elif feature_type == "gene":
                attributes = fields[8]
                gene_name = extract_attribute(attributes, "Name")
                if gene_name:
                    fields[2] = f"gene_{gene_name}"

            adjusted_annotations.append("\t".join(fields))
    
    with open(output_file, 'w') as output:
        for annotation in adjusted_annotations:
            output.write(annotation + "\n")


def fill_unannotated_regions_gff3(annotation_file, dna_length, output_file):
    """Add unannotated region records to a GFF3 annotation file."""

    adjusted_annotations = []
    
    with open(annotation_file, 'r') as ann_file:
        prev_line = None
        for line in ann_file:
            if line.startswith("#") or len(line.strip()) == 0:
                adjusted_annotations.append(line.strip())
                continue
            
            fields = line.strip().split("\t")
            
            region_start = int(fields[3])
            region_end = int(fields[4])

            if prev_line is None and region_start != 1:
                new_fields = fields.copy()
                new_fields[2] = f'UAR_{fields[2]}'
                new_fields[3] = '1'
                new_fields[4] = f'{region_start - 1}'
                adjusted_annotations.append("\t".join(new_fields))
            elif prev_line is not None:
                prev_fields = prev_line.strip().split("\t")
                prev_start = int(prev_fields[3])
                prev_end = int(prev_fields[4])

                if region_start > prev_end + 1:
                    new_fields = fields.copy()
                    new_fields[2] = f'{prev_fields[2]}_UAR_{fields[2]}'
                    new_fields[3] = f'{prev_end + 1}'
                    new_fields[4] = f'{region_start - 1}'
                    adjusted_annotations.append("\t".join(new_fields))

            adjusted_annotations.append("\t".join(fields))
            prev_line = line

    ## last line
    prev_fields = prev_line.strip().split("\t")
    prev_start = int(prev_fields[3])
    prev_end = int(prev_fields[4])
    if dna_length > prev_end + 1:
        new_fields = fields.copy()
        new_fields[2] = f'{prev_fields[2]}_UAR'
        new_fields[3] = f'{prev_end + 1}'
        new_fields[4] = f'{dna_length}'
        adjusted_annotations.append("\t".join(new_fields))

    with open(output_file, 'w') as output:
        for annotation in adjusted_annotations:
            output.write(annotation + "\n")


def extract_attribute(attributes_str, key):
    """
    Extracts the value of a specific attribute from the GFF3 attributes column.
    
    Parameters:
    - attributes_str: The attributes column of the GFF3 line (a semicolon-separated list).
    - key: The key of the attribute to extract (e.g., "product" or "Name").
    
    Returns:
    - The value of the attribute if found, otherwise None.
    """
    attributes = attributes_str.split(";")
    for attribute in attributes:
        if attribute.startswith(f"{key}="):
            return attribute.split("=")[1]
    return None


def add_gaps_before(aligned_seq, pos):
    count = 0
    original_pos = 0
    for base in aligned_seq:
        if base == '-':
            count += 1
        else:
            original_pos += 1
        if original_pos == pos:
            break
    return pos + count


def main():
    parser = argparse.ArgumentParser(description="A tool to align sequences and adjust annotations for gaps.")
    
    parser.add_argument("--reference_fasta", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--sequences_fasta", required=True, help="Path to the sequences FASTA file to align.")
    parser.add_argument("--annotation_file", required=True, help="Path to the annotation GFF3 file.")
    parser.add_argument("--alignment_file", required=True, help="Path to save the aligned sequences FASTA file.")
    parser.add_argument("--clean_alignment_file", required=True, help="Path to save the cleaned alignment FASTA file.")
    parser.add_argument("--output_annotation_file", required=True, help="Path to save the adjusted annotation GFF3 file.")

    args = parser.parse_args()

    run_mafft_alignment_linux(args.reference_fasta, args.sequences_fasta, args.alignment_file)
    remove_gaps_and_ns_from_alignment(args.alignment_file, args.clean_alignment_file)
    adjust_annotation_for_gaps_gff3(args.annotation_file, args.clean_alignment_file, args.output_annotation_file)
    alignment = AlignIO.read(args.clean_alignment_file, "fasta")
    alignment_length = alignment.get_alignment_length()
    fill_unannotated_regions_gff3(args.output_annotation_file, alignment_length, args.output_annotation_file)


if __name__ == "__main__":
    main()

