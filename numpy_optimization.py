from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBParser import PDBParser
from os import listdir
import Bio.PDB.NeighborSearch
import Bio
from Bio.pairwise2 import align
import random
import sys


class CustomModel(Model):
    """Allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity.set_parent(self)
        self.child_list.append(entity)

    def save_to_mmCIF(self, out_name):
        io = Bio.PDB.MMCIFIO()
        io.set_structure(self)
        io.save(out_name+".cif")


class CustomChain(Chain):
    def __init__(self, chainObject):
        self.child_list = chainObject.child_list
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.interactions = []
        self.xtra = chainObject.xtra
        self.level = chainObject.level

    def add_interaction_lst(self, lst):
        self.interactions = lst

    def get_interactions(self, other_chain):
        atomList_1 = list(self.get_atoms())
        atomList_2 = list(other_chain.get_atoms())
        ns = Bio.PDB.NeighborSearch(atomList_2)
        interaction_res_1 = set()
        interaction_res_2 = set()
        for atom in atomList_1:
            interAtoms = ns.search(atom.coord, 5)
            if len(interAtoms) > 0:
                interaction_res_1.add(atom.get_parent().id[1])
                for iatom in interAtoms:
                    interaction_res_2.add(iatom.get_parent().id[1])
        return tuple(sorted(interaction_res_1)), tuple(sorted(interaction_res_2))

    def get_sequence(self):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}
        dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
        rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
        seq = ""
        flag = "prot"
        first_residue_name = self.child_list[0].resname.strip()
        if first_residue_name not in d:
            if "D" in first_residue_name:
                flag = "dna"
            else:
                flag = "rna"
        if flag == "prot":
            for res in self:
                if res.id[0] == " ":
                    seq += d[res.resname]
        elif flag == "dna":
            for res in self:
                if res.id[0] == " ":
                    seq += dna[res.resname.strip()]
        else:
            for res in self:
                if res.id[0] == " ":
                    seq += rna[res.resname.strip()]
        return seq

    def get_common_atoms(self,fixed):
        chain_atoms = sorted(self.get_atoms())
        fixed_atoms = sorted(fixed.get_atoms())
        len_chain = len(chain_atoms)
        len_fixed = len(fixed_atoms)
        if len_chain > len_fixed:
            return chain_atoms[:len_fixed], fixed_atoms
        elif len_fixed > len_chain:
            return chain_atoms, fixed_atoms[:len_chain]
        else:
            return chain_atoms, fixed_atoms

def read_pdbs(directory):
    """Reads pdb Model instances"""
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbmodels = [parser.get_structure("Model_pair", directory + f)[0] for f in listdir(directory)]
    return pdbmodels


def get_new_id(iterator):
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for l in letters:
        if l not in iterator:
            return l
    sys.stderr.write("Too many different chains given. The program can handle modeling maximum 66 different sequences")
    exit(1)


def has_homolgs(target_seq, known_seqs):
    for k_seq in known_seqs:
        alignment = align.globalxx(target_seq, k_seq)[0]
        aln_seq_1 = alignment[0]
        aln_seq_2 = alignment[1]
        al_length = len(alignment[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
        if ident/al_length >= 0.95:
            return k_seq


def unify_ids(pdbmodels):
    """Unifies chain identifiers and returns Custom Model instances"""
    seq_dict = dict()
    for i in range(len(pdbmodels)):
        pdb = pdbmodels[i]
        model = CustomModel(str(i))
        for chain in pdb:
            chain = CustomChain(chain)
            chain.parent = None
            chain_seq = chain.get_sequence()
            if chain_seq not in seq_dict:
                if not seq_dict:
                    new_id = get_new_id(seq_dict.values())
                    seq_dict[chain_seq] = new_id
                    chain.id = new_id
                else:
                    sequences = seq_dict.keys()
                    homolog = has_homolgs(chain_seq, sequences)
                    if homolog:
                        seq_dict[chain_seq] = seq_dict[homolog]
                        chain.id = seq_dict[homolog]
                    else:
                        new_id = get_new_id(seq_dict.values())
                        seq_dict[chain_seq] = new_id
                        chain.id = new_id
            else:
                chain.id = seq_dict[chain_seq]
            model.add(chain)
        pdbmodels[i] = model


def get_interaction_dict(clean_pdbs):
    """Generates interaction dictionary"""
    interaction_dict = dict()
    for pdb in clean_pdbs:
        chain1, chain2 = list(pdb.get_chains())
        inter1_2, inter2_1 = chain1.get_interactions(chain2)
        if inter1_2 != ():
            interaction_dict.setdefault(chain1.id, dict())[inter1_2] = (chain1, chain2, inter2_1)
        if inter2_1 != ():
            interaction_dict.setdefault(chain2.id, dict())[inter2_1] = (chain2, chain1, inter1_2)
    return interaction_dict


def update_interactions_dict(interaction_dict):
    for chain in interaction_dict:
        for interaction_tple in interaction_dict[chain]:
            chain1, chain2, ref_inter = interaction_dict[chain][interaction_tple]
            chain1_filtered_interactions_lst = [x for x in interaction_dict[chain1.id].keys() if x != interaction_tple]
            chain1.add_interaction_lst(chain1_filtered_interactions_lst)
            chain2_filtered_interactions_lst = [x for x in interaction_dict[chain2.id].keys() if x != ref_inter]
            chain2.add_interaction_lst(chain2_filtered_interactions_lst)
            parent = chain1.get_parent()
            parent.child_list = [chain1, chain2]
            interaction_dict[chain][interaction_tple] = chain1, chain2


def has_clashes(move_atoms, model):
    backbone = {"CA", "C1\'"}
    chain_atoms = [atom for atom in move_atoms if atom.id in backbone]
    model_atoms = [atom for atom in model.get_atoms() if atom.id in backbone]
    ns = Bio.PDB.NeighborSearch(model_atoms)
    clashes = 0
    for atom in chain_atoms:
        clashes += bool(ns.search(atom.coord, 2))
    if clashes/len(chain_atoms) >= 0.03:
        return True
    else:
        return False


def get_starting_model(interaction_dict):
    max_len = 0
    interaction_key = None
    for key in interaction_dict:
        l = len(interaction_dict[key])
        if l > max_len:
            max_len = l
            interaction_key = key
    return interaction_dict[interaction_key][next(iter(interaction_dict[interaction_key]))][0].parent.copy()

def build_macrocomplex(directory):
    pdbmodels = read_pdbs(directory)
    unify_ids(pdbmodels)
    interaction_dict = get_interaction_dict(pdbmodels)
    update_interactions_dict(interaction_dict)
    macrocomplex = get_starting_model(interaction_dict)
    run = True
    while run:
        counter = 0
        empty_chains = 0
        for chain in macrocomplex:
            if counter < 300:
                if len(chain.interactions) > 0:
                    random.shuffle(chain.interactions)
                    for inter_tple in chain.interactions:
                        fix, to_move = interaction_dict[chain.id][inter_tple]
                        sup = Bio.PDB.Superimposer()
                        chain_atoms, fix_atoms = chain.get_common_atoms(fix)
                        sup.set_atoms(chain_atoms, fix_atoms)
                        move = to_move.copy()
                        move_atoms = sorted(move.get_atoms())
                        sup.apply(move_atoms)
                        if not has_clashes(move_atoms, macrocomplex):
                            print("Chain " + str(counter) + " added")
                            move.parent = None
                            macrocomplex.add(move)
                            counter += 1
                        else:
                            print("Chain NOT added")
                    chain.interactions = []
                else:
                    print("Empty chain")
                    empty_chains += 1
            else:
                run = False
                break
        if empty_chains == len(macrocomplex):
            run = False
    macrocomplex.save_to_mmCIF(directory.replace("/", ""))


build_macrocomplex("enterovirus/")