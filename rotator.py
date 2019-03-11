from Bio.PDB import *

import numpy

# Parse the structure file

structure = PDBParser().get_structure("DC", "hemo/DC.pdb")[0]

# Iterate through all atoms and rotate by 90 degress

rotation_matrix = rotaxis2m(90, Vector(-9,20,44))

translation_matrix = numpy.array((0, 0, 0), 'f')

for atom in structure.get_atoms():

    atom.transform(rotation_matrix, translation_matrix)

io = PDBIO()

io.set_structure(structure)

io.save("hemo/DC_rotated.pdb")