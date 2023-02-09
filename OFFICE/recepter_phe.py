import sys
from schrodinger import structure
from schrodinger.infra import phase

def main(path):
    with structure.StructureReader(path[0]) as reader:
        for st in reader:
            print(st)            

if __name__ == '__main__':
    main(sys.argv[1:2])
