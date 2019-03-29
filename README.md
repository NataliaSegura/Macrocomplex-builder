# 4SMacro-Builder  
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [4SMacroBuilder](#4smacrobuilder)
- [Software Requirements](#software-requirements)
- [Download and Installation](#download-and-installation)
  - [Installation via PIP](#installation-via-pip)
- [Input Files](#input-files)
- [Tutorial](#tutorial)
  - [Command line arguments](#Command-line-arguments) 
  - [GUI](#gui)
- [Analysis of examples](#Analysis of examples)
  - [Proteosome](#Proteosome)
  - [Enterovirus](#Enterovirus)
  - [Nucleosome](#Nucleosome) 
- [Strong Points](#strong-points)
- [Limitations](#Limitations) 
- [Next steps](#Next-steps)
- [FAQS](#FAQS)
  
<!-- /TOC -->



### 4SMacroBuilder

4SMacroBuilder is a stand-alone python3 program developed by Pau Badia i Monpel, Altaïr C. Hernández and Natàlia Segura Alabart. It builds protein macrocomplexes taking a set of protein-protein, protein-RNA, protein-DNA, RNA - DNA, RNA - RNA, and/or DNA - DNA interactions. This software could serve to study quaternary structures that are difficult to study *in vivo*.

Below is shown how to install and use this program as a standalone command line script (executing the script *MacroB.py*) or with the Graphical User Interface (*tkinter.py*).


### Software Requirements

These are the software and its versions required for the 4SMacroBuilder functionality and execution:

  * [Python 3.6](https://www.python.org/downloads/)
  * [Biopython](http://biopython.org/wiki/Download)
  * argparse
  * os
  * sys 
  * numpy 
  * [Pymol](https://pymol.org/2/)

For the GUI the following ones are also necessary:

  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)


## Download and Installation

You can download our package using Git with the next command. We also recommend creating a directory named "Models":
 
```bash
  $ git clone https://github.com/NataliaSegura/Macrocomplex-builder.git
  $ cd 4SMacroBuilder
  $ mkdir Models
 ```
At this pont, the directory 4SMacroBuilder should contain the following files and directories:

* README.md, README.pdf: the files containing the tutorial and information about our application.
* MBlauncher.py: the command-line script to launch the program.
* Tkinter.py: to launch the GUI.
* MacroB.py: a module requiered by MBlauncher.py where are defined the classes of the program.
* CustomPDB.py: a module required by MacroB.py where are defined the functions of the program.
* Tests: a directory with several examples stored in sub-directories that serve as input to the program.
* Models: an empty folder where the created complexes will be saved.
* Templates: the raw PDB files from which we extracted the example pairwise interactions.
* setup.py script: to install the program in the python side-packages.

Check that all this information has been correctly downloaded and that there is the script called *setup.py*.

In order to be able to use all the scprits provided in 4SMacroBuilder the user has to install the package in the python site-packages.

```bash
   $ python3 setup.py install
```
Be sure to have the dependencies previously stated.


### Installation via PIP

** To be discused **

### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, it can handle up to 95% of identity. 
* The same sequence appearing in different files with the same or different names. 
* Repeated chain interactions are not requiered as inputs (i.e. interaction A-A 10 times is treated as a single PDB). It solves infinite structures (i.e. Microtubul).


## Tutorial

In this section we make a brief explanation of how to use 4SMacroBuilder.

#### Command line Standalone application


```bash
    $ MBlauncher.py -h

    usage: MBlauncher.py [-h] -i, -input DIRECTORY [-o, -output OUTPUT]
                     [-c MAX_CHAINS] [-n NUM_MODELS] [-d] [-v]

    MacrocomplexBuilder is a python program designed to generate macrocomplex
    structures from simple pair inetractions

    optional arguments:
      -h, --help            show this help message and exit
      -i, -input DIRECTORY  Input directory where the pair interaction pdbs are
                            located
      -o, -output OUTPUT    Name of the output file, no extension is needed
      -c MAX_CHAINS         Maximum number of chains that the user wants in the
                            model
      -n NUM_MODELS         Number of models that the program will compute
      -d, --dirty           Generates an output file for each added chain to track
                            how the program builds the complex
      -v, --verbose         Shows what the program is doing
```

### GUI

Another way to use the program is using the the GUI. To do so run the following command:


```bash
$ Tkinter.py
```
For a detailed explanation of how to use the GUI check the *report.pdf*

To get a better understanding of how to run the programme properly, we show different examples that represent different inputs that may be provided. The main aspects that may differ the inputs are: number of different chain interactions and number of atoms of the whole macrocomplex.


## Analysis of examples

Using the template argument. The program can compare between the model it's created and a given template. Using this we can analyse the quality of some reconstructed macrocomplexes.

number of template differences, conta el nombre de diferents cadenes. 
parser.add_argument('-t', dest='template', action="store", type=str, default=None, help="To discriminate against different models, a template can be given to calculate the RMSD")

### Proteosome

1pma is a proteosome from *Thermoplasma acidophilum* (https://www.rcsb.org/structure/1PMA). A proteosome is a protein macrocomplex which degrade proteins.
Throught our template analysis we do not see any differences with respect the number of chains between the one created and the template. Further analysis is done by structural visual comparison using an interactive visualization and molecular structural analysis programm such as UCSF CHIMERA ( image ***XX***). We can't oberve any differences in their structure, the superposition done is totally perfect. 
<p align="center">
  <img src="/proteosome.png" width="450"/>
</p>


### Enterovirus

3j23 is the Enterovirus 71 empty capsids (https://www.rcsb.org/structure/3j23). 
Giving a set of protein-protein interactions, 4SMacroBuilder is able to construct the whole capsid macrocomplex. Comparing the structural composition of the model versus a template, we can observe any differences between both(image ***XX AND YY***).   
<p align="center">
  <img src="/3j23_enterovirus.png" alt="template image" width="450"/>
  <img src="/enterovirus.png" alt="model created" width="450"/>
</p>

### Nucleosome



## Strong Points

* The algorithm is implemented in such a way that the final output is given in a very fast way. This is due to use an interaction check-list while reconstructing the macrocomplex (*see documentation*). 
* The program accepts a wide range of input types: 
  - The input names does not affect to the output (i.e. if all PDB files are named XY.pdb).
  - The input does not need all the interactions in different PDB files (i.e. case of virus capside or microtubul, with more than 150 chain-interaction in the case of the virus capside, and infinite interactions in the microtubul).
  - If the user gives a non-existing or wrong interaction the program ignores it and keeps going.
  - Different ways of executing 4SMAcroBuilder: from command line or with GUI.
  - Possibility to generate different models in a very short time. This allow the user to compare each model and decide which is/are the best models. Different models are scored in an output file.
  - Possibility to model DNA/RNA - protein interactions if given.
  - Possibility to limit the number of chains when executing the program.
  - If the user specifies that wants the macrocomplex with 7 chains but in fact the model has only 4 chains (i.e. Hemoglobin), it will not try to put more just because it was asked. This limited and reduces very much the program performance time.


## Limitations


One of the main limitations dealing with the creation of a macrocomplex is the number of atoms and number of interactions it has. That's why we did a deeper anaylisis of these two factors using the microtubul folder. What is advantatgeous about this macrocomplex is that without any limitations it can go on forever without stopping, more or less like in a cell. But, limiting its parameters, it enables us to analyse our program.

We did a series of test normalizing by number of atoms and interactions. The microtuble has two different chains, with an average of 3347 atoms and 4 interactions by chain. 

As it can be seen in the following graph, the program follows an exponential curve. The more atoms/iterations it has to check, the more time it needs to run.

Due to some aspects of our approach, some "infinite" structures as the microtubule is not modeled as expected. This could be possibly due to a random behavior when superposing different chains. 

We don't apply any optimization method, so there are no option to optimze our model once it's finished. 

## Next Steps

## FAQS
