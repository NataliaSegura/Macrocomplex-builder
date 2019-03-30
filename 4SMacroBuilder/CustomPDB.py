from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
import Bio.PDB.NeighborSearch
import random

class CustomModel(Model):
    """Allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity.set_parent(self)
        self.child_list.append(entity)

    def save_to_mmCIF(self, out_name, verbose):
        io = Bio.PDB.MMCIFIO()
        io.set_structure(self)
        io.save(out_name+".cif")
        if verbose:
            print(out_name+".cif saved")


class CustomChain(Chain):
    """This class adds new behaviours to the biopython's Chain class"""
    protein = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}
    dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
    rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}

    def __init__(self, chainObject):
        self.child_list = chainObject.child_list
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.interactions = []
        self.xtra = chainObject.xtra
        self.level = chainObject.level

    def add_interaction_lst(self, lst):
        """Sets the interactions attribute to a list of interaction tuples"""
        self.interactions = lst

    def get_interactions(self, other_chain):
        """Compares the distance between the atoms of two chains and returns a tuple (chain1, chain2) of
        interaction tuples (2,34,50...)"""
        atom_list_1 = list(self.get_atoms())
        atom_list_2 = list(other_chain.get_atoms())
        ns = Bio.PDB.NeighborSearch(atom_list_2)
        interaction_res_1 = set()
        interaction_res_2 = set()
        for atom in atom_list_1:
            inter_atoms = ns.search(atom.coord, 5)
            if len(inter_atoms) > 0:
                interaction_res_1.add(atom.get_parent().id[1])
                for iatom in inter_atoms:
                    interaction_res_2.add(iatom.get_parent().id[1])
        return tuple(sorted(interaction_res_1)), tuple(sorted(interaction_res_2))

    def get_sequence(self):
        """Returns the chain's sequence, it can be a protein, DNA or RNA sequence"""
        seq = ""
        flag = "prot"
        first_residue_name = self.child_list[0].resname.strip()
        if first_residue_name not in self.protein:
            if "D" in first_residue_name:
                flag = "dna"
            else:
                flag = "rna"
        if flag == "prot":
            for res in self:
                if res.id[0] == " ":
                    seq += self.protein[res.resname]
        elif flag == "dna":
            for res in self:
                if res.id[0] == " ":
                    seq += self.dna[res.resname.strip()]
        else:
            for res in self:
                if res.id[0] == " ":
                    seq += self.rna[res.resname.strip()]
        return seq

    def get_common_atoms(self, other):
        """Compares the list of atoms of two chains and returns an even tuple of atoms"""
        self_atoms = sorted(self.get_atoms())
        other_atoms = sorted(other.get_atoms())
        len_self = len(self_atoms)
        len_other = len(other_atoms)
        if len_self > len_other:
            return self_atoms[:len_other], other_atoms
        elif len_other > len_self:
            return self_atoms, other_atoms[:len_self]
        else:
            return self_atoms, other_atoms