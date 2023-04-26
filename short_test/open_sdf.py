from openbabel import *

obconversion = OBConversion()
obconversion.SetInFormat("sdf")
obmol = OBMol()

notatend = obconversion.ReadFile(obmol,"'/Users/song-inhyeok/Documents/data/all-sdf.sdf'")

print(notatend)