import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_internal_sequence(complete_path, ltr5_path, ltr3_path, output_path):
    # Read the sequences
    complete_seq = SeqIO.read(complete_path, "fasta")
    ltr5_seq = SeqIO.read(ltr5_path, "fasta")
    ltr3_seq = SeqIO.read(ltr3_path, "fasta")

    # Find the start and end positions of the LTRs in the complete sequence
    start_5ltr = str(complete_seq.seq).index(str(ltr5_seq.seq))
    end_5ltr = start_5ltr + len(ltr5_seq)
    start_3ltr = str(complete_seq.seq).rindex(str(ltr3_seq.seq))

    # Extract the internal sequence
    internal_seq = complete_seq.seq[end_5ltr:start_3ltr]

    # Create a new SeqRecord for the internal sequence
    internal_record = SeqRecord(
        Seq(internal_seq),
        id=complete_seq.id + "_internal",
        description="Internal sequence between 5' and 3' LTRs"
    )

    # Write the internal sequence to a new FASTA file
    SeqIO.write(internal_record, output_path, "fasta")
    print(f"Internal sequence extracted and saved to {output_path}")
    print(f"Length of internal sequence: {len(internal_seq)} bp")

def main():
    parser = argparse.ArgumentParser(description="Extract internal sequence between LTRs")
    parser.add_argument("complete", help="Path to complete sequence FASTA file")
    parser.add_argument("ltr5", help="Path to 5' LTR sequence FASTA file")
    parser.add_argument("ltr3", help="Path to 3' LTR sequence FASTA file")
    parser.add_argument("output", help="Path to output FASTA file")
    
    args = parser.parse_args()
    
    extract_internal_sequence(args.complete, args.ltr5, args.ltr3, args.output)

if __name__ == "__main__":
    main()
