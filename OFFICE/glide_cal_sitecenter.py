import sys
from schrodinger import structure
from schrodinger.structutils import analyze

def main(path):
    with structure.StructureReader(path[0]) as reader: 
        for st in reader:
            
            atoms = analyze.evaluate_asl(st, path[1])
            xsum = 0.0
            ysum = 0.0
            zsum = 0.0

            for anum in atoms:
                xsum += st.atom[anum].x
                ysum += st.atom[anum].y
                zsum += st.atom[anum].z

            natoms = float(len(atoms))
            xcent = xsum / natoms
            ycent = ysum / natoms
            zcent = zsum / natoms

            grid_center = [xcent,ycent,zcent] 
            
            print(','.join(str(e) for e in grid_center)) 

if __name__ == '__main__':
    main(sys.argv[1:3])
