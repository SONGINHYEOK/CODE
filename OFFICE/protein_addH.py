import sys
from schrodinger import structure
from schrodinger.structutils import build

def main(path):
    with structure.StructureReader(path[0]) as reader:
        for st in reader:
            build.add_hydrogens(st)      
            with structure.StructureWriter(path[1]) as writer:
                writer.append(st)
    
	
if __name__ == '__main__':
    main(sys.argv[1:3])
