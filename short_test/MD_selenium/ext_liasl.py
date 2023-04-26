from schrodinger import structure
from schrodinger.structutils import analyze

st = structure.StructureReader.read("/Users/song-inhyeok/CODING/short_test/align_1V4S-out.cms")
ligands = analyze.find_ligands(st)


ligand = ligands[0]
print(ligand.ligand_asl)




