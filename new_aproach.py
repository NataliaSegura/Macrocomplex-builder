from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBParser import PDBParser
from os import listdir
from Bio.PDB.Polypeptide import PPBuilder
import Bio.PDB.NeighborSearch
import copy

class CustomModel(Model):
    """Allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity_id = entity.get_id()
        entity.set_parent(self)
        self.child_list.append(entity)
        self.child_dict[entity_id] = entity


class CustomChain(Chain):
    def __init__(self, chainObject, lst):
        self.child_dict = chainObject.child_dict
        self.child_list = chainObject.child_list
        self.full_id = chainObject.full_id
        self.level = chainObject.level
        self.xtra = chainObject.xtra
        self._id = chainObject._id
        self.parent = chainObject.parent
        ppb = PPBuilder()
        self.sequence = ppb.build_peptides(chainObject)[0].get_sequence()._data
        self.interactions = lst


    def add_interaction_lst(self, lst):
        pass


def read_pdbs(directory):
    """Reads pdb Model instances"""
    parser = PDBParser(PERMISSIVE=1)
    pdbfiles = [directory + f for f in listdir(directory)]
    pdbpairs = []
    for pdbf in pdbfiles:
        pdbpairs.append(parser.get_structure("Model_pair", pdbf)[0])
    return pdbpairs

def get_new_id(iterator):
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for l in letters:
        if l not in iterator:
            return l

def get_seq_dict(pdbmodels):
    """Unifies chain identifiers and returns Custom Model instances"""
    seq_dict = dict()
    ppb = PPBuilder()
    models = []
    i = 0
    for pdb in pdbmodels:
        model = CustomModel(str(i))
        i += 1
        chain1, chain2 = list(pdb.get_chains())
        chain1.detach_parent()
        chain2.detach_parent()
        chain1_seq = ppb.build_peptides(chain1)[0].get_sequence()._data
        if chain1_seq not in seq_dict:
            new_id = get_new_id(seq_dict.values())
            seq_dict[chain1_seq] = new_id
            chain1.id = new_id
        else:
            chain1.id = seq_dict[chain1_seq]
        chain2_seq = ppb.build_peptides(chain2)[0].get_sequence()._data
        if chain2_seq not in seq_dict:
            new_id = get_new_id(seq_dict.values())
            seq_dict[chain2_seq] = new_id
            chain2.id = new_id
        else:
            chain2.id = seq_dict[chain2_seq]
        model.add(chain1)
        model.add(chain2)
        models.append(model)
    return seq_dict, models

def get_interactions(chain1, chain2):
    atomList_1 = list(chain1.get_atoms())
    atomList_2 = list(chain2.get_atoms())
    ns = Bio.PDB.NeighborSearch(atomList_2)
    interaction_res_1 = set()
    interaction_res_2 = set()
    for atom in atomList_1:
        interAtoms = ns.search(atom.coord, 3.5)
        if len(interAtoms) > 0:
            interaction_res_1.add(atom.get_parent().id[1])
            for iatom in interAtoms:
                interaction_res_2.add(iatom.get_parent().id[1])
    return tuple(sorted(interaction_res_1)), tuple(sorted(interaction_res_2))

def get_interaction_dict(clean_pdbs):
    """Generates interaction dictionary"""
    interaction_dict = dict()
    for pdb in clean_pdbs:
        chain1, chain2 = list(pdb.get_chains())
        inter1_2, inter2_1 = get_interactions(chain1, chain2)
        interaction_dict.setdefault(chain1.id, dict())[inter1_2] = (chain1, chain2)
        interaction_dict.setdefault(chain2.id, dict())[inter2_1] = (chain2, chain1)
    return interaction_dict

def update_interactions_dict(interaction_dict):
    """Updates each interaction's chains with a CustomChain class"""
    for chain in interaction_dict:
        for interaction_tple in interaction_dict[chain]:
            chain1, chain2 = interaction_dict[chain][interaction_tple]
            chain1 = CustomChain(chain1, list(interaction_dict[chain1.id].keys()))
            chain2 = CustomChain(chain2, list(interaction_dict[chain2.id].keys()))
            parent = chain1.get_parent()
            parent.child_list = [chain1, chain2]
            interaction_dict[chain][interaction_tple] = chain1, chain2

def has_clashes(chain, model):
    chain_atoms = list(chain.get_atoms())
    model_atoms = list(model.get_atoms())
    ns = Bio.PDB.NeighborSearch(model_atoms)
    clashes = 0
    for atom in chain_atoms:
        clashes += len(ns.search(atom.coord, 1.5))

    if clashes/len(chain_atoms) >= 0.05 :
        return True
    else:
        return False


directory = "hemo/"
raw_pdbmodels = read_pdbs(directory)
seq_dict, clean_pdbs = get_seq_dict(raw_pdbmodels)
interaction_dict = get_interaction_dict(clean_pdbs)
update_interactions_dict(interaction_dict)

run = True
new_pdb = interaction_dict["A"][(104, 106, 107, 108, 109, 181, 184, 398)][0].get_parent()

limit = 0
while run:
    counter = 0
    for chain in new_pdb:
        if limit < 8:
            if len(chain.interactions) > 0:
                for inter_tple in chain.interactions:
                    fix, move = interaction_dict[chain.id][inter_tple]
                    sup = Bio.PDB.Superimposer()
                    sup.set_atoms(list(chain.get_atoms()), list(fix.get_atoms()))
                    sup.apply(list(move.get_atoms()))
                    if not has_clashes(move, new_pdb):
                        new_pdb.add(copy.deepcopy(move))
                    x = 0
                chain.interactions = []
                counter = 1
                limit += 1
        else:
            run = False
        if counter != 0:
            run = False

io = Bio.PDB.PDBIO()
io.set_structure(new_pdb)
io.save(directory + 'macrocomplex.pdb')
x = 0


