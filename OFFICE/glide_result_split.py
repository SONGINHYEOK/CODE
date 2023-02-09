import sys
from schrodinger import structure

def main(path):
    with structure.StructureReader(path[0]) as reader:
        for i,st in enumerate(reader):
            if i !=0:
                with structure.StructureWriter('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/docking/'+str(i)+'_'+path[1]) as writer:
                    writer.append(st)
                    
if __name__ == '__main__':
    main(sys.argv[1:3])
