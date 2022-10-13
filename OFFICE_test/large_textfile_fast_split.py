from itertools import islice

lines_per_file = 5000
with open("/Users/song-inhyeok/Documents/ZINC20_cr_project/final_filtering_ZINC20.smi") as file:
    i = 1
    while True:
        try:
            checker = next(file)
        except StopIteration:
            break
        with open(f"/Users/song-inhyeok/Documents/ZINC20_cr_project/final_ZINC20/ZINC20_{i}.smi", 'w') as out_file:
            out_file.write(checker)
            for line in islice(file, lines_per_file-1):
                out_file.write(line)
        i += 1
