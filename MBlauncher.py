from MacroB import build_macrocomplex
from argparse import ArgumentParser
import sys

parser = ArgumentParser(description='MacrocomplexBuilder is a python program designed to generate macrocomplex structures from simple pair inetractions')
stoich = parser.add_mutually_exclusive_group()
parser.add_argument("-i", dest="directory", action="store", type=str, help="Input directory where the pair interaction pdbs are located", required=True)
parser.add_argument("-o", dest="output", action="store", type=str, default='macrocomplex.pdb', help="Name of the output file, no extension is needed")
parser.add_argument('-c', dest='max_chains', action="store", type=int, default=300, help="Maximum number of chains that the user wants in the model")
parser.add_argument('-n', dest='num_models', action="store", type=int, default=1, help="Number of models that the program will compute")
parser.add_argument('-d', '--dirty', dest='dirty', action="store_true", default=False, help="Generates an output file for each added chain to track how the program builds the complex")
parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", default=False, help="Shows what the program is doing")
stoich.add_argument('-t', dest='template', action="store", type=str, default=False, help="To discriminate against different models, a template can be given to calculate the RMSD")
stoich.add_argument('-s', dest='stoich_string', action='store', type=str, default=False, help="The user can also give the desired stechiometry in this format: A:6,B:11,C:2 ...")
options = parser.parse_args()
build_macrocomplex(options.directory, options.output, options.max_chains, options.num_models, dirty=options.dirty,
                       verbose=options.verbose, template=options.template, stoich_string=options.stoich_string)


x = 0