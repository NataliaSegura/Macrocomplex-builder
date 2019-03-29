# 4SMacro-Builder  
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Index]
- [4SMacroBuilder](#4smacrobuilder)
- [Software Requirements](#software-requirements)
- [Download and Installation](#download-and-installation)
  - [Installation via PIP](#installation-via-pip)
- [Input Files](#input-files)
- [Tutorial](#tutorial)
  - [Command line arguments](#Command-line-arguments) 
  - [GUI](#gui)
- [Examples](#Examples)
  - [Example 1](#Example-1)
  - [Example 2](#Example-2)
  - [Example 3](#Example-3) 
  - [Example 4](#Example-3) 
  - [Example 5](#Example-3)  
- [Strong Points](#strong-points)
- [Limitations](#Limitations) 
- [Next steps](#Next-steps)
- [FAQS](#FAQS)
  
<!-- /TOC -->



### 4SMacroBuilder

4SMacroBuilder is an stand-alone python3 program developed by Pau Badia i Monpel, Altaïr C. Hernández and Natàlia Segura Alabart that builds protein macrocomplexes taking a set of protein-protein interactions. This software is also able to construct macrocomplexes with DNA or RNA interactions, and could serve to study quaternary structures that are difficult to study *in vivo*.

Below is shown how to install and use this program as a standalone command line script(executing the script *MacroB.py*) or with the Graphical User Interface (*tkinter.py*).



### Software Requirements

To run 4SMacroBuilder with all its functionalities some software versions are required, those we used to test the program.

  * [Python 3.6](https://www.python.org/downloads/)
  * [Biopython](http://biopython.org/wiki/Download)
  * argparse
  * os
  * sys 
  * numpy 

For the GUI the following ones are also necessary:

  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)
  * [Pymol](https://pymol.org/2/)


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
* MacroB.py: a module requiered by MBlauncher.py where are defined the classes of the program.
* CustomPDB.py: a module required by MacroB.py where are defined the functions of the program.
* Tests: a directory with several examples stored in sub-directories that serve as input to the program.
* Models: an empty folder where the created complexes will be saved.
* Templates: the raw PDB files from which we extracted the example pairwise interactions.
* setup.py script

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

* The same sequence appearing in different PDB files has not to be identical, we can handle 95% of identity. 
* The same sequence appearing in different files can have different names. 
* Repeated chain interactions are not requiered as inputs (i.e. interaction A-A 10 times is treated as a single PDB). It solves        infinite structures (i.e. Microtubul).


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

Another way to use the program is using the the GUI. To do so just run the following command:


```bash

$ Tkinter.py

```
For a detailed explanation of how to use the GUI check the *report.pdf*

To get a better understanding of how to run the programme properly, we show different examples that represent different inputs that may be provided. The main aspects that may differ the inputs are: number of different chain interactions and number of atoms of the whole macrocomplex.


## Examples

### Example 1

### Example 2

### Example 3

### Example 4

### Example 5


## Strong Points

## Limitations

## Next Steps

## FAQS
