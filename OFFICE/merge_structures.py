from __future__ import print_function
import sys
from schrodinger import structure

if len(sys.argv)<3 or '-h' in sys.argv:
    print("")
    print("usage: $SCHRODINGER/run merge_structures.py <infile> [<infile2>...] <outfile>")
    print("")
    print("       Merge all structures from one or multiple files into a")
    print("       single structure written to <outfile>.")
    print("")
    sys.exit()

infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

print("")
print("Merging structures from:")
for infile in infiles:
    print("    %s" % infile)
print("")
print("Output file: %s" % outfile)

nst = 0
outst = None
for infile in infiles:
   for st in structure.StructureReader(infile):
       nst += 1
       if outst==None:
           outst = st
       else:
           outst = outst.merge(st,copy_props=True)

outst.write(outfile)
print("")
print("Done. Merged %s structures." % nst)
