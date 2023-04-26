from schrodinger import structure
from schrodinger.structutils import analyze

lig_info = {}

import os 

file_list = os.listdir("/Users/song-inhyeok/Documents/PDB_prep/ori")
count = 0
for i in file_list:
    if i != ".DS_Store":
        count +=1
        with structure.StructureReader(f'/Users/song-inhyeok/Documents/PDB_prep/ori/{i}') as reader:
            for st in reader:    
                ligands = analyze.find_ligands(st)
                res, *rest = st.residue
                if len(ligands)==0:
                    with open("./sch_not_ligand_list.text", "a") as f:
                        f.write(i.split('.')[0]+ '\n')
                        f.close()
                else:
                    continue
        print(count)