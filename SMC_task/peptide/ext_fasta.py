import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def process_id_seq_pair(start, step, fasta_lines):
    id_lines = fasta_lines[start::step]
    seq_lines = fasta_lines[start + 1::step]

    # Remove '>' from the ID lines
    id_lines = [line.strip('>\n') for line in id_lines]

    return list(zip(id_lines, seq_lines))

def save_id_seq_pairs(input_fasta, output_csv, num_processes=2):
    with open(input_fasta, 'r') as fasta_file:
        fasta_lines = fasta_file.readlines()

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # Process odd lines with one set of workers
        future_odd = executor.submit(process_id_seq_pair, 0, 2, fasta_lines)
        
        # Process even lines with another set of workers
        future_even = executor.submit(process_id_seq_pair, 1, 2, fasta_lines)

        # Get results from futures
        id_seq_pairs_odd = future_odd.result()
        id_seq_pairs_even = future_even.result()

    # Combine results and create a DataFrame excluding lines containing "Delete"
    id_seq_pairs_all = id_seq_pairs_odd + id_seq_pairs_even
    id_seq_pairs_filtered = [(id, seq.rstrip()) for id, seq in id_seq_pairs_all if "Delete" not in seq]

    # Create a DataFrame from the filtered pairs
    df = pd.DataFrame(id_seq_pairs_filtered, columns=['ID', 'seq'])

    # Save DataFrame to CSV
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    input_fasta_file = "./APD_Hs_all.fasta"
    output_csv_file = "./human_all_peptide.csv"

    save_id_seq_pairs(input_fasta_file, output_csv_file)

    print(f"ID and seq pairs (excluding lines with 'Delete') saved to {output_csv_file}")
