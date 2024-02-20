import pandas as pd
import os

def csv_to_fasta(csv_filename, output_folder):
    # Read the CSV file using pandas
    df = pd.read_csv(csv_filename)

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # Assuming the column containing peptide sequences is named 'peptide sequence'
        peptide_sequence = row['peptide sequence'].strip()

        # Create a unique filename for each sequence
        fasta_filename = os.path.join(output_folder, f'sequence_{index + 1}.fasta')

        # Write the FASTA header and sequence to the file
        with open(fasta_filename, 'w') as fasta_file:
            fasta_file.write(f'>{peptide_sequence}\n{peptide_sequence}\n')

# Example usage
csv_filename = '/Users/song-inhyeog/CODEING/CODE/test.csv'  # Replace with your CSV file path
output_folder = '/Users/song-inhyeog/CODEING/CODE/test_result'  # Replace with the desired output folder path

csv_to_fasta(csv_filename, output_folder)
