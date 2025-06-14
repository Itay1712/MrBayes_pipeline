from Bio import AlignIO
from collections import defaultdict
import argparse

def create_nexus_from_alignment_and_annotation(alignment_file, annotation_file, output_nexus):
    """
    Creates a NEXUS file with partition information based on the alignment and GFF3 annotations.
    
    Parameters:
    - alignment_file: Path to the multiple sequence alignment (in FASTA format).
    - annotation_file: Path to the GFF3 annotation file.
    - output_nexus: Path to save the generated NEXUS file.
    """
    
    # Step 1: Parse the alignment file (in FASTA format) using Biopython
    alignment = AlignIO.read(alignment_file, "fasta")
    sequences = {record.id: str(record.seq) for record in alignment}
    ntax = len(sequences)   # Number of taxa (sequences)
    nchar = len(alignment[0])  # Number of characters (alignment length)
    
    # Step 2: Parse the GFF3 file to extract regions and calculate region lengths
    region_lengths = defaultdict(int)
    regions = []
    
    with open(annotation_file, 'r') as gff:
        for line in gff:
            if line.startswith("#") or len(line.strip()) == 0:
                continue
            fields = line.strip().split("\t")
            region_name = fields[2]  # Feature type (e.g., gene, exon)
            region_start = int(fields[3])
            region_end = int(fields[4])
            region_length = region_end - region_start + 1
            region_lengths[region_name] += region_length
            regions.append((region_name, region_start, region_end))
    
    # Step 3: Generate NEXUS file
    with open(output_nexus, "w") as nexus_file:
        # Write NEXUS header
        nexus_file.write("#NEXUS\n\n")
        nexus_file.write("BEGIN DATA;\n")
        nexus_file.write(f"    DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        nexus_file.write("    FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
        nexus_file.write("    MATRIX\n")
        
        # Write concatenated sequences
        for sample, sequence in sequences.items():
            nexus_file.write(f"{sample}    {sequence}\n")
        
        nexus_file.write("    ;\n")
        nexus_file.write("END;\n\n")
        
        # Write partition information using the GFF3 regions
        nexus_file.write("BEGIN SETS;\n")
        for region_name, region_start, region_end in regions:
            nexus_file.write(f"    CHARSET {region_name} = {region_start}-{region_end};\n")
        nexus_file.write("END;\n\n")
        # Write MrBayes block for partitioned analysis
        nexus_file.write("BEGIN MRBAYES;\n")

        nexus_file.write("\n")
        nexus_file.write("    set autoclose=yes;\n")
        nexus_file.write("    lset nst=6 rates=gamma;\n")
        
        partition_mtdna = ", ".join([f"{region_name}" for region_name in region_lengths.keys()])
        nexus_file.write(f"    partition mtDNA = {len(region_lengths)}: {partition_mtdna};\n")
        nexus_file.write("    set partition=mtDNA;\n")
        nexus_file.write("    unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);\n")
        nexus_file.write("    prset ratepr=variable;\n")
        nexus_file.write("    mcmcp ngen=1000000 samplefreq=1000 nchains=4 savebrlens=yes;\n")
        nexus_file.write("    mcmc;\n")
        nexus_file.write("    sump burnin=250000;\n")
        nexus_file.write("    sumt burnin=250000;\n")
        nexus_file.write("END;\n")

    print(f"NEXUS file has been created: {output_nexus}")


def main():
    parser = argparse.ArgumentParser(description="A tool to create a Nexus file from sequence alignment and annotation.")

    parser.add_argument("--alignment_file", required=True, help="Path to the aligned sequences FASTA file.")
    parser.add_argument("--annotation_file", required=True, help="Path to the annotation GFF3 file.")
    parser.add_argument("--output_nexus", required=True, help="Path to save the output Nexus file.")

    args = parser.parse_args()

    create_nexus_from_alignment_and_annotation(args.alignment_file, args.annotation_file, args.output_nexus)

if __name__ == "__main__":
    main()


