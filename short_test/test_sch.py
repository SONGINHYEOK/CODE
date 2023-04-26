from schrodinger import structure
from schrodinger.structutils import analyze
import sys            

with structure.StructureReader("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/1V4S_prep.maegz") as reader: 
        for st in reader:
            
            atoms = analyze.evaluate_asl(st, 'res.num 501')
            xsum = 0.0
            ysum = 0.0
            zsum = 0.0

            for anum in atoms:
                xsum += st.atom[anum].x
                ysum += st.atom[anum].y
                zsum += st.atom[anum].z

            natoms = float(len(atoms))
            if xsum !=0:
                xcent = xsum / natoms
            else:
                xcent=0.0
            if ysum !=0:
                ycent = ysum / natoms
            else:
                ycent = 0.0
            if zsum !=0:
                zcent = zsum / natoms
            else:
                zcent = 0.0
            grid_center = [xcent,ycent,zcent] 
            
            print(','.join(str(e) for e in grid_center)) 
        