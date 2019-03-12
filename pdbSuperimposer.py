import Bio.PDB
import numpy
from Bio.PDB.PDBParser import PDBParser
from os import listdir
from os.path import isfile, join
from Bio.PDB.Chain import Chain
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Model import Model
from Bio.PDB.Entity import Entity
import string

parser = PDBParser(PERMISSIVE=1)

directory = "phosphate/"
pdbfiles = [directory + f for f in listdir(directory) if isfile(join(directory, f))]


# alt_model = parser.get_structure("CB", "CB.pdb")[0]["A"]
class CustomModel(Model):
    def add(self, entity):
        """Add a child to the Entity."""
        entity_id = entity.get_id()
        entity.set_parent(self)
        self.child_list.append(entity)
        self.child_dict[entity_id] = entity



class Interaction(object):
    def __init__(self, tple, structure, parentId, model):
        self.tple = tple
        self.structure = structure
        self.parentId = parentId
        self.model = model

    def __hash__(self):
        return hash(tuple(self.structure))

    def __eq__(self, other):
        return self.structure == other.structure


class ChainStructure(Chain):
    def __init__(self, chainObject):
        chainObject.detach_parent()
        self.child_dict = chainObject.child_dict
        self.child_list = chainObject.child_list
        self.full_id = chainObject.full_id
        self.level = chainObject.level
        self.xtra = chainObject.xtra
        self._id = chainObject._id
        self.parent = chainObject.parent
        ppb = PPBuilder()
        self.sequence = ppb.build_peptides(chainObject)[0].get_sequence()._data
        self.interactions = set()

    def __hash__(self):
        return self.sequence.__hash__()

    def __eq__(self, other):
        return self.sequence == other

    def get_interactions(self, chainstructe):
        atomList_1 = self.get_atoms()
        atomList_2 = list(chainstructe.get_atoms())
        ns = Bio.PDB.NeighborSearch(atomList_2)
        interaction_atoms_1 = []
        interaction_atoms_2 = set()
        for atom in atomList_1:
            interAtoms = ns.search(atom.coord, 7.5)
            if len(interAtoms) > 0:
                interaction_atoms_1.append(atom)
                for a in interAtoms:
                    interaction_atoms_2.add(a)
        pair = CustomModel(self.id+chainstructe.id)
        pair.add(self)
        pair.add(chainstructe)
        return Interaction((self.id, chainstructe.id), sorted(interaction_atoms_1), self.id, pair), \
               Interaction((chainstructe.id, self.id), sorted(interaction_atoms_2), chainstructe.id, pair)

    def set_inter(self, interaction):
        self.interactions.add(interaction)




class Macrocomplex(object):
    def __init__(self):
        self.chains = {}
        self.idIndex = 0

    def add_pair(self, model):
        pair = []
        for chain in model:
            chain = ChainStructure(chain)
            id = self.chain_exists(chain)
            if not id:
                new_id = self.make_id()
                chain.id = new_id
                self.chains[new_id] = chain
            else:
                chain.id = id
            pair.append(chain)
        chain1, chain2 = pair
        ichain1, ichain2 = chain1.get_interactions(chain2)

        id = self.chain_exists(chain1)
        if id in self.chains:
            self.chains[id].set_inter(ichain1)
        else:
            if chain1.id == chain2.id:
                chain1.set_inter(ichain1)
                chain1.set_inter(ichain2)
                self.chains[chain1.idChain] = chain1
            else:
                chain1.set_inter(ichain1)
                chain2.set_inter(ichain2)
                self.chains[chain1.idChain] = chain1
                self.chains[chain2.idChain] = chain2
        id = self.chain_exists(chain2)
        if id in self.chains:
            self.chains[id].set_inter(ichain2)
        else:
            if chain1.id == chain2.id:
                chain1.set_inter(ichain1)
                chain1.set_inter(ichain2)
                self.chains[chain1.idChain] = chain1
            else:
                chain1.set_inter(ichain1)
                chain2.set_inter(ichain2)
                self.chains[chain1.idChain] = chain1
                self.chains[chain2.idChain] = chain2

    def chain_exists(self, chain):
        for key, seq in self.chains.items():
            if seq == chain.sequence:
                return key
        return False

    def make_id(self):
        id = string.ascii_uppercase[self.idIndex]
        self.idIndex += 1
        return id

def get_sequence(chain):
    ppb = PPBuilder()
    return ppb.build_peptides(chain)[0].get_sequence()._data


def read_pdb_file(pdbf, m):
    model = parser.get_structure("Model_pair", pdbf)[0]
    m.add_pair(model)
    x= 0


def read_pdb_dir(file_list):
    m = Macrocomplex()
    for f in file_list:
        read_pdb_file(f, m)
    return m


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


m = read_pdb_dir(pdbfiles)

def recursive(m, parent_k, pairs, macrocomp, used_keys):
    for child_k in dict[parent_k]:
        if macrocomp == []:
            used_keys.append(parent_k)
            macrocomp = dict[parent_k][child_k]
            recursive(dict, child_k, pairs, macrocomp, used_keys)
        else:
            if child_k not in used_keys:
                parent_chain = dict[parent_k][child_k][parent_k]
                parent_atoms = sorted(parent_chain.get_atoms())
                # Superimpose#
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
io.save(directory + 'macrocomplex.pdb')


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
