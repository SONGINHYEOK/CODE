import pandas as pd
import difflib
import argparse
from multiprocessing import Pool

def calculate_similarity(args):
    input_string, target_string, similarity_threshold = args
    # Get the ratio of similarity between two strings
    similarity_ratio = difflib.SequenceMatcher(None, input_string, target_string).ratio()
    
    # Check if similarity meets the threshold
    if similarity_ratio >= similarity_threshold:
        return target_string, similarity_ratio

    # Return None if similarity is below the threshold
    return None

def search_and_save_results(input_string, data_frame, column_name, output_csv_path, num_processes, similarity_threshold):
    search_results = []
    similarity_ratios = []
    indices = []

    # Create a list of arguments for the calculate_similarity function
    args_list = [(input_string, str(target_string), similarity_threshold) for target_string in data_frame[column_name]]

    # Use multiprocessing to parallelize the calculations
    with Pool(num_processes) as pool:
        results = pool.map(calculate_similarity, args_list)

    # Extract the results
    for i, result in enumerate(results):
        if result is not None:
            search_results.append(result[0])
            similarity_ratios.append(result[1])
            indices.append(i)

    # Create a new DataFrame with search results, similarity ratios, and indices
    results_df = pd.DataFrame({
        'Search_Result': search_results,
        'Similarity_Ratio': similarity_ratios,
        'Index': indices
    })

    # Save the results to a CSV file
    results_df.to_csv(output_csv_path, index=False)
    print(f"Search results saved to {output_csv_path}")

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description='Search and calculate similarity with CSV data.')
    parser.add_argument('input_string', help='Input string to search for')
    parser.add_argument('--csv_path', default='/Users/song-inhyeog/CODEING/CODE/SMC_task/peptide/peptide_DB_complex.csv', help='Path to CSV file')
    parser.add_argument('--column_name', default='Peptide Sequence', help='Name of the column to search in CSV')
    parser.add_argument('--output_csv', default='./search_results.csv', help='Path to save the search results CSV')
    parser.add_argument('--num_processes', type=int, default=2, help='Number of processes for multiprocessing')
    parser.add_argument('--similarity_threshold', type=float, default=0.6, help='Similarity threshold')

    args = parser.parse_args()

    # Read CSV file
    data_frame = pd.read_csv(args.csv_path)

    # Perform search and save results to a CSV file
    search_and_save_results(args.input_string, data_frame, args.column_name, args.output_csv, args.num_processes, args.similarity_threshold)
