from Bio.PDB.PDBParser import PDBParser
from os import listdir
import Bio.PDB.NeighborSearch
from Bio.pairwise2 import align
import random
import sys
from CustomPDB import CustomModel, CustomChain


def read_pdbs(directory):
    """Reads the input directory and generates pdb models"""
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbmodels = [parser.get_structure("Model_pair", directory + f)[0] for f in listdir(directory)]
    return pdbmodels


def get_new_id(iterator):
    """Returns a new pdb id that is not in the given iterator"""
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for l in letters:
        if l not in iterator:
            return l
    sys.stderr.write("Too many different chains given. The program can handle modeling maximum 66 different sequences")
    exit(1)


def has_homolgs(target_seq, known_seqs):
    """Checks if a given sequence is an homolog of any of the known sequences"""
    for k_seq in known_seqs:
        alignment = align.globalxx(target_seq, k_seq)[0]
        aln_seq_1 = alignment[0]
        aln_seq_2 = alignment[1]
        al_length = len(alignment[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
        if ident/al_length >= 0.95:
            return k_seq


def unify_ids(pdbmodels):
    """Unifies chain identifiers and updates the pdb models to CustomModel class"""
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
    """Updates the interactions attribute of each chain inside the dictionary"""
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
    """Compares the atoms backbone atoms of the moving chain with the backbone atoms of the model"""
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


def get_starting_model(interaction_dict, verbose):
    """Returns as a starting model the CustomModel inside any of the interactions of the chain with more interactions"""
    max_len = 0
    interaction_key = None
    for key in interaction_dict:
        length = len(interaction_dict[key])
        if length > max_len:
            max_len = length
            interaction_key = key
    inter_tple = next(iter(interaction_dict[interaction_key]))
    return interaction_dict[interaction_key][inter_tple][0].parent.copy()


def build_macrocomplex(directory, output, max_chains, num_models, dirty, verbose):
    """Main funtion"""
    pdbmodels = read_pdbs(directory)
    unify_ids(pdbmodels)
    interaction_dict = get_interaction_dict(pdbmodels)
    update_interactions_dict(interaction_dict)
    for i in range(num_models):
        print("Building macrocomplex number "+str(i)+":")
        macrocomplex = get_starting_model(interaction_dict, verbose)
        run = True
        counter = 2
        empty_chains = 0
        while run:
            for chain in macrocomplex:
                if counter < max_chains:
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
                                if verbose:
                                    print("Chain " + str(counter) + " added: interaction "+chain.id+": " +
                                          str(inter_tple[0:2])+" ... "+ str(inter_tple[-2:]) + " to "+move.id)
                                move.parent = None
                                macrocomplex.add(move)
                                counter += 1
                                if dirty:
                                    macrocomplex.save_to_mmCIF(output+"_tmp_"+str(counter))
                            elif verbose:
                                print("Chain NOT added: interaction " + chain.id + ": " +
                                      str(inter_tple[0:2]) + " ... " + str(inter_tple[-2:]) + " to " + move.id)
                        chain.interactions = []
                    else:
                        print("Empty chain")
                        empty_chains += 1
                else:
                    run = False
                    break
            if empty_chains >= len(macrocomplex):
                run = False
        print("Macrocomplex number "+str(i)+" finished")
        macrocomplex.save_to_mmCIF(output+"_"+str(i), verbose)
