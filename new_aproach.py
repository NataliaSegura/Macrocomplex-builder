from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBParser import PDBParser
from os import listdir
from Bio.PDB.Polypeptide import PPBuilder
import Bio.PDB.NeighborSearch
import Bio
import copy

class CustomModel(Model):
    """Allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity_id = entity.get_id()
        entity.set_parent(self)
        self.child_list.append(entity)
        #self.child_dict[entity_id] = entity


class CustomChain(Chain):
    def __init__(self, chainObject, lst):
        self.child_dict = chainObject.child_dict
        self.child_list = chainObject.child_list
        self.full_id = chainObject.full_id
        self.level = chainObject.level
        self.xtra = chainObject.xtra
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.sequence = get_sequence(self)
        self.interactions = lst


    def add_interaction_lst(self, lst):
        pass


def read_pdbs(directory):
    """Reads pdb Model instances"""
    parser = PDBParser(PERMISSIVE=1)
    pdbfiles = [directory + f for f in listdir(directory)]
    pdbpairs = []
    for pdbf in pdbfiles:
        model = parser.get_structure("Model_pair", pdbf)[0]
        for chain in model:
            seq = get_sequence(chain)
            chain_res_list = []
            for res in chain:
                if res.id[0] == ' ':
                    chain_res_list.append(res)
            chain.child_list = chain_res_list
        pdbpairs.append(model)
    return pdbpairs

def get_new_id(iterator):
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for l in letters:
        if l not in iterator:
            return l

def get_seq_dict(pdbmodels):
    """Unifies chain identifiers and returns Custom Model instances"""
    seq_dict = dict()
    models = []
    i = 0
    for pdb in pdbmodels:
        model = CustomModel(str(i))
        i += 1
        chain1, chain2 = list(pdb.get_chains())
        chain1.detach_parent()
        chain2.detach_parent()
        chain1_seq = get_sequence(chain1)
        if chain1_seq not in seq_dict:
            new_id = get_new_id(seq_dict.values())
            seq_dict[chain1_seq] = new_id
            chain1.id = new_id
        else:
            chain1.id = seq_dict[chain1_seq]
        chain2_seq = get_sequence(chain2)
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
        if inter1_2 != ():
            interaction_dict.setdefault(chain1.id, dict())[inter1_2] = (chain1, chain2)
        if inter2_1 != ():
            interaction_dict.setdefault(chain2.id, dict())[inter2_1] = (chain2, chain1)
    return interaction_dict

def get_sequence(chain):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    seq = ""
    for res in chain:
        seq += d[res.resname]
    return seq


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

def has_clashes(move_atoms, model):
    backbone = {"CA"}
    chain_atoms = []
    for atom in move_atoms:
        if atom.id in backbone:
            chain_atoms.append(atom)
    model_atoms = []
    for atom in model.get_atoms():
        if atom.id in backbone:
            model_atoms.append(atom)
    ns = Bio.PDB.NeighborSearch(model_atoms)
    clashes = 0
    for atom in chain_atoms:
        clashes += len(ns.search(atom.coord, 1.5))

    if clashes/len(chain_atoms) >= 0.015:
        return True
    else:
        return False


directory = "mosaic_virus/"
raw_pdbmodels = read_pdbs(directory)
seq_dict, clean_pdbs = get_seq_dict(raw_pdbmodels)
interaction_dict = get_interaction_dict(clean_pdbs)
update_interactions_dict(interaction_dict)


#new_pdb = interaction_dict["A"][(104, 106, 107, 108, 109, 181, 184)][0].get_parent()
#new_pdb = interaction_dict["A"][(34, 37, 40, 43, 99, 101, 105, 146)][0].get_parent()
#new_pdb = interaction_dict["A"][(122, 124, 132, 133, 135, 139, 140, 141, 157, 164, 165)][0].get_parent()
#new_pdb = interaction_dict["A"][(79, 128, 129, 130, 146, 150, 151, 152, 153, 154, 155, 156, 158, 162, 199, 200, 201, 203, 204, 205, 206)][0].get_parent()

new_pdb = interaction_dict["A"][next(iter(interaction_dict["A"]))][0].get_parent().copy()


run = True
limit = 0
while run:
    counter = 0
    for chain in new_pdb:
        if counter < 9999999:
            if len(chain.interactions) > 0:
                for inter_tple in chain.interactions:
                    fix, to_move = interaction_dict[chain.id][inter_tple]
                    sup = Bio.PDB.Superimposer()
                    sup.set_atoms(sorted(chain.get_atoms()), sorted(fix.get_atoms()))
                    move = to_move.copy()
                    move_atoms = sorted(move.get_atoms())
                    sup.apply(move_atoms)
                    if not has_clashes(move_atoms, new_pdb):
                        print("Chain " + str(counter) + " added")
                        move.parent = None
                        new_pdb.add(move)
                        io = Bio.PDB.PDBIO()
                        io.set_structure(new_pdb)
                        io.save('micro/' + str(counter) + '.pdb')
                        counter += 1
                    else:
                        print("Chain NOT added")
                    x = 0
                    limit += 1
                chain.interactions = []
            else:
                print("Empty chain")
        else:
            run = False
            break
    if counter != 0:
        print("Done")
        run = False
"""
for chain in interaction_dict:
    for interaction in interaction_dict[chain]:
        model = CustomModel(chain+str(interaction))
        chain1, chain2 = interaction_dict[chain][interaction]
        model.add(chain1)
        model.add(chain2)
        io = Bio.PDB.PDBIO()
        io.set_structure(model)
        io.save(chain+str(interaction) + '.pdb')"""
"""
io = Bio.PDB.PDBIO()
io.set_structure(new_pdb)
io.save('macrocomplex_'+ directory.replace("/", "_") +'.pdb')"""
final_file = CustomModel("Macrocomplex")
for chain in new_pdb:
    chain.detach_parent()
    chain.id = get_new_id([x.id for x in final_file.get_chains()])
    final_file.add(chain)


io = Bio.PDB.PDBIO()
io.set_structure(final_file)
io.save('macrocomplex_'+ directory.replace("/", "_") +'.pdb')
x = 0

io = Bio.PDB.MMCIFIO()
io.set_structure(final_file)
io.save(directory.replace("/", "_")+".cif")
