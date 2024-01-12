import streamlit as st
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from multiprocessing import Pool

def calculate_sequence_similarity(args):
    input_sequence, target_sequence, similarity_threshold = args
    # Convert the input sequences to Bio.Seq objects
    peptide_seq1 = Seq(input_sequence)
    peptide_seq2 = Seq(target_sequence)

    # Use the Needleman-Wunsch algorithm for global sequence alignment
    alignments = pairwise2.align.globalxx(peptide_seq1, peptide_seq2, one_alignment_only=True, score_only=True)

    # Calculate the sequence similarity as a percentage
    sequence_similarity = (alignments / max(len(peptide_seq1), len(peptide_seq2)))

    # Check if similarity meets the threshold
    if sequence_similarity >= similarity_threshold:
        return sequence_similarity

    # Return None if similarity is below the threshold
    return None

def search_and_save_results(input_sequence, data_frame, column_name, similarity_threshold):
    similarity_scores = []

    # Create a list of arguments for the calculate_sequence_similarity function
    args_list = [(input_sequence, str(target_sequence), similarity_threshold) for target_sequence in data_frame[column_name]]

    # Use multiprocessing to parallelize the calculations
    with Pool() as pool:
        # Use pool.map with the calculate_sequence_similarity function
        results = pool.map(calculate_sequence_similarity, args_list)

    # Create a new DataFrame with calculated similarity scores
    results_df = pd.DataFrame({
        'Sequence_Similarity': results
    })

    # Merge the calculated similarity column with the original data
    merged_df = pd.concat([data_frame, results_df], axis=1)

    return merged_df

# Streamlit UI
def main():
    st.title("Peptide Sequence Similarity Search App")

    # User input
    input_sequence = st.text_input("Enter the search sequence:")

    # Read CSV file
    csv_file_path = st.file_uploader("Upload CSV file", type=["csv"])
    if csv_file_path is not None:
        data_frame = pd.read_csv(csv_file_path)
        st.success("CSV file uploaded successfully!")
    else:
        st.warning("Please upload a CSV file.")
        st.stop()

    # Choose column to search
    column_name = st.selectbox("Select the column to search in:", data_frame.columns, index=data_frame.columns.get_loc("Peptide Sequence"))

    # Set similarity threshold
    similarity_threshold = st.slider("Set similarity threshold", 0.0, 1.0, 0.8, step=0.01)

    # Perform search and display results
    if st.button("Search"):
        results_df = search_and_save_results(input_sequence, data_frame, column_name, similarity_threshold)

        # Display only rows meeting the similarity threshold
        st.write("Search results:")
        st.write(results_df[results_df['Sequence_Similarity'] >= similarity_threshold])

        # Save all rows to CSV
        st.markdown("### Save Results")
        save_path = st.text_input("Enter the path to save the results CSV:")
        if st.button("Save Results") and not save_path.isspace():
            results_df.to_csv(save_path, index=False)
            st.success(f"Results saved to {save_path}")

if __name__ == "__main__":
    main()
