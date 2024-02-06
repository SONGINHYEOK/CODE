from admet_ai import ADMETModel

model = ADMETModel()
preds = model.predict(smiles="O(c1ccc(cc1)CCOC)CC(O)CNC(C)C")

import csv

def save_dict_to_csv(file_path, data):
    # Extracting the keys to use them as column headers
    headers = data.keys()

    # Writing the dictionary to a CSV file
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        
        # Writing the header
        writer.writeheader()

        # Writing the data
        writer.writerow(data)

# Example usage:

save_dict_to_csv('output.csv', preds)