with open('/Users/song-inhyeok/Downloads/Chemspace_Freedom_Space_Set/Chemspace_Freedom_Space_SMILES.smiles') as myfile:
    total_lines = sum(1 for line in myfile)

print(total_lines)