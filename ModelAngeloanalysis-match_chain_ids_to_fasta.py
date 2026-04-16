# Made by Luca Troman October 2024
# Made to take the cif file output from model angelo and match this to the fasta file you input to model angelo.
# Run with python thisscript.py --model path/to/your/model.cif --fasta path/to/your/sequences.fasta --output path/to/your/output.csv

import argparse
import csv
from Bio import SeqIO
from Bio.PDB import MMCIFParser

# Dictionary to convert three-letter amino acid codes to single-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def convert_three_to_one(three_letter_seq):
    """Convert a three-letter amino acid sequence to a single-letter sequence."""
    single_letter_seq = ''
    for i in range(0, len(three_letter_seq), 3):
        res = three_letter_seq[i:i+3]
        single_letter_seq += three_to_one.get(res, '-')  # Append a gap for unrecognized sequences
    return single_letter_seq

def is_subsequence(part, full):
    """Check if `part` can be a subsequence of `full` allowing gaps."""
    part = part.replace('-', '')  # Remove gaps from part
    full = full.replace('-', '')  # Remove gaps from full
    m, n = len(part), len(full)
    
    # Use a two-pointer technique to allow for gaps
    i, j = 0, 0
    while i < m and j < n:
        if part[i] == full[j]:  # Characters match
            i += 1
        j += 1  # Move in the full sequence always

    return i == m  # If we reached the end of 'part', it's a match

def get_simplified_name(header, mapping):
    """Get the simplified protein name based on the FASTA header."""
    return mapping.get(header, 'Not found')

def main(cif_file, fasta_file, output_csv):
    # Load the CIF file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('model', cif_file)

    # Extract chain IDs and their amino acid sequences
    chain_data = {}
    for model in structure:
        for chain in model:
            amino_acids = ''.join([residue.get_resname() for residue in chain])
            chain_data[chain.id] = amino_acids

    # Load the FASTA file into a dictionary
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    # Define the mapping of FASTA headers to simplified protein names
    protein_name_mapping = {
        'tr|Q73M00|Q73M00_TREDE': 'FlaA1',
        'tr|Q73MU8|Q73MU8_TREDE': 'FlaA2',
        'tr|Q73K73|Q73K73_TREDE': 'HEAT_domain',
        'tr|Q73MN1|Q73MN1_TREDE': 'FlaB1',
        'tr|Q73NZ6|Q73NZ6_TREDE': 'FlaB2',
        'tr|Q73MN3|Q73MN3_TREDE': 'FlaB3',
        'tr|Q73MU9|Q73MU9_TREDE': 'FlaA3',
    }

    # Prepare to write to the output CSV
    with open(output_csv, mode='w', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(['Chain ID', 'Amino Acids (3-letter)', 'Amino Acids (1-letter)', 'Amino Acid Count', 'FASTA Header', 'Simplified Protein Name'])  # Header row

        for chain_id, amino_acids in chain_data.items():
            # Convert to single-letter
            single_letter_seq = convert_three_to_one(amino_acids)

            # Calculate amino acid count
            amino_acid_count = len(single_letter_seq.replace('-', ''))

            # Initially mark the FASTA header and simplified name as not found
            fasta_header = 'Not found'
            simplified_name = 'Not found'

            # Check if amino acids match any FASTA sequence
            for header, sequence in fasta_dict.items():
                if is_subsequence(single_letter_seq, sequence):
                    fasta_header = header
                    simplified_name = get_simplified_name(header, protein_name_mapping)
                    break  # Stop after finding a match
            
            # Write the data to the CSV
            writer.writerow([chain_id, amino_acids, single_letter_seq, amino_acid_count, fasta_header, simplified_name])

    print(f"Mapping and matching saved to {output_csv}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Map chains from a CIF file to sequences in a FASTA file and convert amino acids.")
    parser.add_argument('--model', required=True, help='Path to the CIF file (model).')
    parser.add_argument('--fasta', required=True, help='Path to the FASTA file (sequences).')
    parser.add_argument('--output', default='chain_sequence_mapping.csv', help='Output CSV file name (default: chain_sequence_mapping.csv).')

    # Parse arguments
    args = parser.parse_args()

    # Run main function with provided arguments
    main(args.model, args.fasta, args.output)

