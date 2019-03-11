import Bio.PDB
import numpy
from Bio.PDB.PDBParser import PDBParser
from os import listdir
from os.path import isfile, join

parser = PDBParser(PERMISSIVE=1)


pdbfiles = ["phosphate/" + f for f in listdir("phosphate/") if isfile(join("phosphate/", f))]
# alt_model = parser.get_structure("CB", "CB.pdb")[0]["A"]

def read_pdb(file_list):
    pairs = []
    for f in file_list:
        pairs.append(parser.get_structure(f.split(".")[0], f)[0])
    return pairs

def get_atom_interactions(chain1, chain2):
    chain1_atoms = list(chain1.get_atoms())
    ns = Bio.PDB.NeighborSearch(chain1_atoms)
    interaction_res_chain1 = set()
    interaction_res_chain2 = set()
    for atom in chain2.get_atoms():
        interAtoms = ns.search(atom.coord, 7.5)
        if len(interAtoms) > 0:
            interaction_res_chain2.add(atom)
            for intatom in interAtoms:
                interaction_res_chain1.add(intatom)
    interactions = sorted(interaction_res_chain1.union(interaction_res_chain2))
    return interactions

pairs = read_pdb(pdbfiles)
interaction = {}
for model in pairs:
    chains = []
    for chain in model:
        chains.append(chain)
    inter = get_atom_interactions(chains[0], chains[1])
    interaction.setdefault(chains[0].id, {}).setdefault(chains[1].id, model)
    interaction.setdefault(chains[1].id, {}).setdefault(chains[0].id, model)


# Get root tree
top_key = ""
c = 0
for key in interaction.keys():
    k_number = len(interaction[key])
    if k_number > c:
        c = k_number
        top_key = key

def recursive(dict, parent_k, pairs, macrocomp, used_keys):
    for child_k in dict[parent_k]:
        if macrocomp == []:
            used_keys.append(parent_k)
            macrocomp = dict[parent_k][child_k]
            recursive(dict, child_k, pairs, macrocomp, used_keys)
        else:
            if child_k not in used_keys:
                parent_chain = dict[parent_k][child_k][parent_k]
                parent_atoms = sorted(parent_chain.get_atoms())
                #Superimpose#
                sup = Bio.PDB.Superimposer()
                sup.set_atoms(sorted(macrocomp[parent_k].get_atoms()), parent_atoms)
                child_chain = dict[parent_k][child_k][child_k]
                sup.apply(child_chain)
                macrocomp.add(child_chain)
                #############
                used_keys.append(parent_k)
                recursive(dict, child_k, pairs, macrocomp, used_keys)
    return macrocomp

macrocomp = []
used_keys = []
macrocomp = recursive(interaction, top_key, pairs, macrocomp, used_keys)

io = Bio.PDB.PDBIO()
io.set_structure(macrocomp)
io.save('phosphate/macrocomplex.pdb')





class InteractionSelect(Bio.PDB.Select):
    def __init__(self, intres):
        self.intres = intres

    def accept_residue(self, residue):
        for atom in residue:
            atomcoords = tuple(atom.coord)
            if atomcoords in self.intres:
                return 1
            else:
                return 0




"""
sup = Bio.PDB.Superimposer()
chainC = parser.get_structure("CA", "hemo/CA.pdb")[0]["C"]
chainA = parser.get_structure("CA", "hemo/CA.pdb")[0]["A"]

CB = parser.get_structure("CD", "hemo/DC_rotated.pdb")[0]

interCA = get_atom_interactions(chainC, chainA)[0]
chainC = CB["C"]
move_atoms = []
for atom in chainC.get_atoms():
    if atom in interCA:
        move_atoms.append(atom)

sup = Bio.PDB.Superimposer()
sup.set_atoms(interCA, move_atoms)
sup.apply(CB)
io = Bio.PDB.PDBIO()
io.set_structure(CB)
io.save('hemo/transformedCD.pdb')
"""
x = 0
