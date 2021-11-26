from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import numpy
import xpdb

'''parser = PDBParser()
structure = parser.get_structure("estructura", "1a1d.pdb")'''
sloppyparser = PDBParser(
    PERMISSIVE=True, structure_builder=xpdb.SloppyStructureBuilder()
)
structure = sloppyparser.get_structure("estructura", "1a1d.pdb")
sloppyio = xpdb.SloppyPDBIO()
##print(sloppyio.get_structure())
##Bio.PDB.Selection.uniqueify(items)

'''for model in structure:
    for chain in model:
        for residue in chain:
            print(residue)'''
for atom in structure.get_atoms():
    print(atom)

for residue in structure.get_residues():
    print(residue)




