import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def filter_sequence_length(row):
    # Function to check if the length of the sequence is 10 or less
    return len(str(row['Peptide Sequence'])) == 10

def parallel_filter(df_chunk):
    # Apply the filter_sequence_length function to each row
    mask = df_chunk.apply(filter_sequence_length, axis=1)
    
    # Apply the mask to the DataFrame to get the filtered results
    filtered_df_chunk = df_chunk[mask]

    return filtered_df_chunk

def main():
    # Replace 'your_data.csv' with the actual filename/path.
    input_file_path = '/Users/song-inhyeog/Downloads/human_all_peptide.csv'
    output_file_path = 'filtered_data.csv'

    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_file_path)

    # Split the DataFrame into chunks based on the number of available CPU cores
    num_cores = 4  # You can adjust this based on your machine's capabilities
    chunk_size = len(df) // num_cores
    chunks = [df.iloc[i:i + chunk_size] for i in range(0, len(df), chunk_size)]

    # Process each chunk in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor() as executor:
        filtered_dfs = list(executor.map(parallel_filter, chunks))

    # Concatenate the results from different chunks
    result_df = pd.concat(filtered_dfs, ignore_index=True)

    # Save the filtered DataFrame to a new CSV file
    result_df.to_csv(output_file_path, index=False)

    print(f"Filtered data saved to {output_file_path}")

if __name__ == "__main__":
    main()
