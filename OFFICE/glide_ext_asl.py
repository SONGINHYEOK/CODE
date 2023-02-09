from schrodinger import structure
from schrodinger.structutils import analyze
import sys
lig_info = {}

def main(path):
    with structure.StructureReader(path[0]) as reader:
        for st in reader:    
            ligands = analyze.find_ligands(st)
            for ligand in ligands:
                name = ligand.pdbres.strip()
                lig_info[name]=ligand.ligand_asl
    print(lig_info[path[1]])


if __name__ == '__main__':
    main(sys.argv[1:3])
