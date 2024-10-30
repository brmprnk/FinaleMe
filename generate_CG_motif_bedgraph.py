import gzip
from Bio import SeqIO
from tqdm import tqdm

def find_cpg_motifs(fasta_file, output_file):
    # Open the gzipped FASTA file
    with gzip.open(fasta_file, "rt") as handle:
        # Parse the FASTA file
        fasta_sequences = SeqIO.parse(handle, "fasta")

        # Open output file for writing BEDGRAPH data
        with open(output_file, "w") as out_bedgraph:
            # Iterate over each sequence record
            for fasta in tqdm(fasta_sequences, desc="Processing sequences"):
                chr_name = fasta.id
                sequence = str(fasta.seq).upper()

                # Optionally uncomment the following line to process only chr22
                # if chr_name != "chr22":
                #     continue

                # Search for CpG motifs
                for i in range(len(sequence) - 1):
                    if sequence[i:i+2] == "CG":
                        # BEDGRAPH format: chrom start end value
                        start = i
                        end = i + 2
                        value = 1  # You can set any value here, usually 1 for presence
                        out_bedgraph.write(f"{chr_name}\t{start}\t{end}\t{value}\n")

if __name__ == "__main__":
    # Input hg19 reference genome in FASTA format (gzipped)
    fasta_file = "hg19.fa.gz"
    # Output bedgraph file
    output_file = "CG_motif_chr22.bedgraph"
    
    find_cpg_motifs(fasta_file, output_file)
