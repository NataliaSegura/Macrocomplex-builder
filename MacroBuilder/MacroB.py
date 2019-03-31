from Bio.PDB.PDBParser import PDBParser
from os import listdir
import Bio.PDB.NeighborSearch
from Bio.pairwise2 import align
import random
import sys
from CustomPDB import CustomModel, CustomChain


def read_pdbs(directory, verbose=False):
    """Reads the input directory and generates pdb models"""
    if verbose:
        print("Reading pdb input files from %s" % directory)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbmodels = [parser.get_structure("Model_pair", directory + f)[0] for f in listdir(directory)]
    if verbose:
        print("Pdb objects stored")
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


def unify_ids(pdbmodels, verbose=False):
    """Unifies chain identifiers and updates the pdb models to CustomModel class"""
    seq_dict = dict()
    if verbose:
        print("Unified Ids")
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
    if verbose:
        print("Ids unified")
    return seq_dict


def get_interaction_dict(clean_pdbs, verbose=False):
    """Generates interaction dictionary"""
    interaction_dict = dict()
    if verbose:
        print("Generating interaction dictionary...")
    for pdb in clean_pdbs:
        chain1, chain2 = list(pdb.get_chains())
        inter1_2, inter2_1 = chain1.get_interactions(chain2)
        if inter1_2 != ():
            interaction_dict.setdefault(chain1.id, dict())[inter1_2] = (chain1, chain2, inter2_1)
        if inter2_1 != ():
            interaction_dict.setdefault(chain2.id, dict())[inter2_1] = (chain2, chain1, inter1_2)
    if verbose:
        print("Interaction dictionary generated")
    return interaction_dict


def update_interactions_dict(interaction_dict, verbose=False):
    """Updates the interactions attribute of each chain inside the dictionary"""
    if verbose:
        print("Updating interaction dictionary...")
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
    if verbose:
        print("Interaction dictionary updated")


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


def get_starting_model(interaction_dict, verbose=False):
    """Returns as a starting model the CustomModel inside any of the interactions of the chain with more interactions"""
    if verbose:
        print("Selecting best interaction from where to start modeling...")
    max_len = 0
    interaction_key = None
    for key in interaction_dict:
        length = len(interaction_dict[key])
        if length > max_len:
            max_len = length
            interaction_key = key
    inter_tple = next(iter(interaction_dict[interaction_key]))
    if verbose:
        print("First two chains added")
    return interaction_dict[interaction_key][inter_tple][0].parent.copy()


def generate_model_profile(model):
    """Generates a dictionary with the id chain as key and the number of repetitions of this chain as values"""
    profile = {}
    for chain in model:
        profile.setdefault(chain.id, 0)
        profile[chain.id] += 1
    return profile


def compare_profiles(model_profile, template_profile):
    """Compares a given model profile with the template's, returns the number of differences"""
    score = 0
    for key, value in model_profile.items():
        if key in template_profile:
            if template_profile[key] > value:
                score += template_profile[key] - value
            else:
                score += value - template_profile[key]
        else:
            score += value
    return score


def revaluate_models(model_lst, template, seq_dict, verbose=False):
    if verbose:
        print("Re-evaluating models according to given template...")
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    template = parser.get_structure("template", template)[0]
    template_profile = {}
    for chain in template:
        chain = CustomChain(chain)
        chain_sequence = chain.get_sequence()
        if chain_sequence in seq_dict:
            template_profile.setdefault(seq_dict[chain_sequence], 0)
            template_profile[seq_dict[chain_sequence]] += 1
    models_profiles = []
    for model in model_lst:
        models_profiles.append((compare_profiles(generate_model_profile(model), template_profile), model))
    models_profiles.sort(key=lambda tup: tup[0])
    if verbose:
        print("Models sorted according to given template")
    return models_profiles


def build_macrocomplex(directory, output, max_chains, num_models, template, dirty, verbose):
    """Main funtion"""
    print("Program is running, please wait...")
    in_pdbmodels = read_pdbs(directory, verbose)
    seq_dict = unify_ids(in_pdbmodels, verbose)
    interaction_dict = get_interaction_dict(in_pdbmodels, verbose)
    update_interactions_dict(interaction_dict, verbose)
    out_models = []
    for i in range(1, num_models+1):
        print("Macrocomplex "+str(i)+" ...")
        macrocomplex = get_starting_model(interaction_dict, verbose).copy()
        macrocomplex.id = "Model_"+str(i)
        run = True
        num_of_chains = 2
        num_empty_chains = 0
        while run:
            for chain in macrocomplex:
                if num_of_chains < max_chains:
                    if chain.interactions:
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
                                    print("Chain " + str(num_of_chains) + " added: interaction "+chain.id+": " +
                                          str(inter_tple[0])+" ... "+ str(inter_tple[-1]) + " to "+move.id)
                                move.parent = None
                                macrocomplex.add(move)
                                num_of_chains += 1
                                if dirty:
                                    macrocomplex.save_to_mmCIF(output+"_tmp_"+str(num_of_chains))
                            elif verbose:
                                print("Chain NOT added: interaction " + chain.id + ": " +
                                      str(inter_tple[:1]) + " ... " + str(inter_tple[-1]) + " to " + move.id)
                        chain.interactions = False
                    else:
                        if verbose:
                            print("Chain "+ chain.id + " empty")
                        num_empty_chains += 1
                else:
                    run = False
                    break
            if num_empty_chains >= len(macrocomplex):
                run = False
        print("Macrocomplex "+str(i)+" finished")
        out_models.append(macrocomplex)
    if template:
        ordered_models = revaluate_models(out_models, template, seq_dict, verbose)
        i = 1
        print("Saving models and generating report.log ...")
        with open("report.log", "w") as f:
            f.write("Model\tNum_diff\n")
            for diff, model in ordered_models:
                model.id = "Model_"+str(i)
                model.save_to_mmCIF(output + "_" + str(i), verbose)
                f.write(output + "_" + str(i)+"\t"+str(diff)+"\n")
                i += 1
    else:
        i = 1
        print("Saving models...")
        for model in out_models:
            model.save_to_mmCIF(output + "_" + str(i), verbose)
            i += 1
    print("Done")