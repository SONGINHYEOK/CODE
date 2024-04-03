import pandas as pd
import itertools

df = pd.read_csv("./pp.csv", index_col=None)

start_list = df['0'].tolist()

#amino_acid = "ARNDCQEGHILKMFPSTWYV"


#    print(st)
 
import itertools

amino_acids = "ARNDCQEGHILKMFPSTWYV"


seq_list = []

# Iterate over each index in st
for st in start_list:
    for i in range(len(st)):
        # Iterate over each character in amino_acids
        for amino_acid in amino_acids:
            # Replace the character at index i with the current amino_acid
            modified_sequence = st[:i] + amino_acid + st[i+1:]
            # Print the modified sequence
            seq_list.append(modified_sequence)

unique_list = list(set(seq_list))

total_df = pd.DataFrame()
total_df['sequence'] = unique_list
total_df.to_csv("./PD_derived_total_seq.csv", index=False)
