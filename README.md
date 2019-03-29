# Macrocomplex-Builder  
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [4SMacroBuilder](#sbi-python-project)
  - [Index](#Index)
  - [Introduction](#Introduction)
  - [Objective](#Objective)
- [Requirements](#Requirements)
- [Installation](#general-information)
- [Approach](#apporach)
  - [Algorithm](#algorithm)
  - [Strong Points](#strong-points)
  - [Limitations](#Limitations) 
- [Tutorial](#tutorial)
  - [Command-line arguments](#Command-line-arguments) 
  - [Input files](#input-files)
  - [Example 1](#Example-1)
  - [Example 2](#Example-2)
  - [Example 3](#Example-3)  
- [GUI - UNDER DEVELOPMENT](#gui-under-development)
  - [Launching the app](#Launching-the-app)
- [Next steps](#Next-steps)
- [FAQ](#FAQS)
  
<!-- /TOC -->



### Introduction

Proteins are the most versatil macromolecules in living systems, as are in charge on multitude of specific chemical trasnformations, which provide the cell with usable energy and the molecules needed to form its structure and mainatin the *homeostasis*. Proteins also recieve signals from outside the cell, starting intracellular signal transductions and regulating the gene expression in different stress situations. Sometimes, proteins monomers interact between other monomers in order to form protein macrocomplexes, known as the quaternary structure of a protein. Some examples are hemoglobin, the ATP syntase, the RNA polymerase or the ribosome. 

Nevertheless, understanding how proteins interact with others in the assembly process is not an easy task. For that reason, different research groups have developed methods that predict how this interaction may occur. In the Protein Data Bank (PDB) are stored a large set of proteins with known structure, after a process of x-ray crystallography or Nuclear Magnetic Ressonance, allowing us to study the molecular space and possible allosteric interactions. 

It is important to stand out that there are some proteins that are difficult to crystallize due to its molecular conformation or dimensionallity, such as virus capsides. In other cases, it could be found a protein - DNA/RNA interaction, and it could be an interesting feature to study in process of transcriptional or translational processes (it would be the case of ribosomes, for instance). 


**Objective**

The main scope of this project is to reconstruct protein macrocomplexes from individual protein pairwise interactions. In order to do so, we developed a stand-alone application that reads a set of protein - protein interactions in PDB format and reconstruct different multi-subunit complexes. 


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
  * [Pymol](https://www.pymol.org/install.html?)


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
* seup.py script

Check that all this information has been correctly downloaded and that there is the script called *setup.py*.

In order to be able to use all the scprits provided in 4SMacroBuilder the user has to install the package in the python site-packages.

```bash
   $ python3 setup.py install
```
Be sure to have the dependencies previously stated.


### Installation via PIP

** To be discused **

## Approach

### Algorithm implementation


The aim of this project is to build a protein macrocomplex (quaternary structure) with a single previous knowledge: a set of know protein pairwise interacions. Consider a protein assembly with uncertainties on the shape and/or the position of its constituting protein subunits. Know let's suppose that we know that there are the following protein - protein interactions: 

  * A-B
  * A-C
  * B-D

  *Here A,B,C and D are chains or subunits of the macrocomplexe and may differ in sequence and structure*. 

We could start taking one of these pair interactions as a template (i.e. A-B), and then superpose the rest of the interactions by protein superposition. In protein macrocomplexes there are several chains that interact with more than one chain, allowing the rest of interactions to be done. Assuming that at least one chain of the template has another interaction (in this example A-C), we can align the simillar chains in order to move the new pair interaction to the template, imposing the condition that the same chains are structurally superposed. Therefore, we would obtain a resulting structure of three chains (C-A-B). The next step would follow the same idea, but this time the tample would be the output of the first iteration (chains C-A-B), so we could keep superposing the simillar chains and moving the new pair interaction to the resulting template in previous iteration. Finally we would obtain the final resulting macrocomplex.

**images**

In order to carry this out we should know the order in which the program would have to superpose these pair of interactions, in order to avoid clashes between chains or even to prevent the program to superimpose the same chain more than one time. Also, we should know how many iterations the program would have to achive in order to make the final structure. 

This could be solved in different ways. For instance, starting with one pair of interactions as template one could check all the possible interactions in each iteration an see which candidate would satisfy the problem statement(exhaustive search algorithm). Although this approach would be simple to implement, it would have a computational cost proportional to the number of candidate solutions, which would tend to grow in an exponential way.  

We approached this limitation basing on the interacting residues of each subunit. If we imagine the whole macrocomplex structure as a lego puzzle we could realize that each chain has some residues that are interacting with at least another chain (hydrophobic residues) and the rest of the residues that are exposed to the solvent environment (hydrophilic residues). If we base on that premise, we could store those interacting residues for each chain, as well as the corresponding chain those residues are interacting with. This is possible as we have this informations in the PDB files, so that would be the first task to do. In that way we could force the program to check in each iteration/superimposition whether those residues are interacting or not. At the same time we would have to consider as feasible complexes those that has no clashes when superimposed, which means that the backbone of the model that is been superimopsed is not interacting with the rest of the complex already joined (this means a threshold distance of 2 Aº).   

Then, our program would start to structurally superimpose structures with at least two identical subunits (those that share a pairwise sequence identity >= 95%), but this time for each model the program knows how many interacting sites are in each protein, and even with which specific chain has to interact on those sites. 

Let's explain this with an example: 

Imagine that we have the same pairwise interactions as before explained:

  * A-B
  * A-C
  * B-D

And know, when we read the PDB files, for each chain we store the interacting residues that take place by looking for those residues that are in a distance no larger that 5 Aº. Besides, we store in the same set which is the model of the corresponding interaction, so that we know which is the model to superimpose later. 

For instance:

               A -> [(1,2,3),(8,9,10)]       -> model A-B
                 -> [(15,16),(22,23,24,25)]  -> model A-C

               B -> [(5,6,9),(13,15,17)]     -> model A-B
                 -> [(24,25,28),(33,34)]     -> model B-D

               C -> [(5,6,9),(20,22)]        -> model A-C

               D -> [(1,2,4,5),(13,15,17)]   -> model B-D








### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, we can handle 95% of identity. 
* The same sequence appearing in different files can have different names. 










## Tutorial

The following section explains the user how to use our application to make the most out of it.

### Command-line arguments

* -i --input: directory containing the PDB files that will serve as input to the program.                        
* -o --output: directory where the macrocomplexes will be saved as PDB files
* -v --verbose: whether or not the user wants to keep track of what the programme is doing.                       
* -st --stoich: the stoichiometry of the multi-subunit complexes.                     
* -opt --optimization: if the user wants to optimize the resulting complexes.        
* -cn --contact_num: number of amino acids that are allowed to protrude after each superimposition (default=5).
* -f --files_num: maximum number of files (complexes) the user wants to obtain.                   

To get a better understanding of how to run the programme properly, we show 3 different scenarios that represent the heterogeneity of inputs that may be provided.

#### Example 1

The first example corresponds to the phosphate dehydratase (pdb id: 2F1D). As we can see in the splited file we provide, (2f1d_split.pdb), this is an Homo 8-mer (stoichiometry: A8). To achieve this complex we can run:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A8" -f10 -v
```
Where example1 is the directory containing all input files, complexes is the directory where the output files will be saved, "A8" is the stoichiometry (i.e. 8 equal chains), -f10 means that the programme will stop after 10 complexes, and -v means that the standard error will be printed. In this case, we can see how multiple complexes are output. This is because at the recurssion level the function stopped (k=8), the programme could still superimpose more structures without finding clashes. Interestingly, we see how the complex9.pdb is the same as the 2f1d_split (inside raw_pdbs), except for one chain.

As more chains can be added, we can rerun the programme with a higher stoichiometry: A24. This time though, it will take longer to execute, so we will limit the max number of files to 1. Before doing so, we must move all the previous complexes to a different file or change the outdir:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A24" -f1 -v
```
As in this example we only have one unique chain, if we try to get a hetero k-mer, we will get an error:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A2B2" -f1 -v
```

#### Example 2

The second example was kindly provided by Prof. Baldo Oliva. In this case, we have two unique chains, so we can build hetero k-mers:

```python
python3 multicomplex.py -i example2 -o complexes -st "A6B6" -f1 -v 
```
As we can see, it completes the recursion pretty fast. However, if we try to add two more chains the resulting complexes will clash, so it will take too long and do not give any result:

```python
python3 multicomplex.py -i example2 -o complexes -st "A7B7" -f1 -v 
```
Thus, if we are not interested in any particular stoichiometry, it is a good practice to start with few recursions, check the output, and increase it or not accordingly. Finally, in the case were we have more than one unique chain, we can still try to find Homo k-mers

```python
python3 multicomplex.py -i example2 -o complexes -st "A6" -f1 -v 
```
#### Example 2 with optimization

This example shows how to run the program with optimization:

```python
python3 multicomplex.py -i example2 -o complexes -st "A6B6" -f1 -v -opt
```
As we can see, the program will create several additional files on the folder from which we are calling the program execution:
* complexX.rsr : this file contains the restraints used to optimize the structure.
* complexX.ini : this file contains the initial MODELLER model.
* complexX.D00000001 and complexX.D9999XXXX.pdb : those files contain the progress of optimization. 
* complexX_optimized.pdb : This file contains the final optimized structure. 

### Example 3

The third example corresponds to the 20S proteasome (pdb id: 1G65). This represents an extreme case, where every chain is different. The proteasome is a Hetero 28-mer (stoichiometry: A2B2C2D2E2F2G2H2I2J2K2L2M2N2). For our purposes, we splited in two and obtained a stoichiometry of: A1B1C1D1E1F1G1H1I1J1K1L1M1N1 (1g65_split, inside raw_pdbs). We will input increasing stoichiometries and assess the scalability of the programme accordingly:

```python
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1" -f4 -v 
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1E1F1G1" -f4 -v 
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1E1F1G1H1I1J1K1L1M1N1" -f4 -v 

```

The last example gives exactly half 20S proteasome. However, it took extremely long to execute (~15min). Thus, if one wants to obtain a complex with a great deal of unique chains, it is a wise idea to let the programme run and go have a coffee ;).


## GUI - UNDER DEVELOPMENT

We are currently developing a GUI interface to run the program without command line. Here I will show the work done and the next steps. 

### Launching the app

To do so, you have to run the following command on the terminal:

```bash
python3 gui2.py
```
The next window will appear on your computer. 

<img src="GUIpics/main_window.png" width="500" height="500">

We can see a the window divided into two frames. On the left frame we will have a list with the files on the specified directory. Then there is a Clean button to restore all the options on the right and the Run program! button to launch the request. If we take a look on the right frame we can see the different options the user can specify. We can chech the optimize en verbose checkboxes, enter the maximum desired files, browse to choose the output directory and finally an Entry to specify the stoichiometry of the desired macrocomplex. 

Finally looking at the upper part of the app there is the menu, with two main parts: "Folder/File Opts" and "Help". 

| *Folder/File Opts* | *Help Menu* |
| :---: | :---: |
| <img src="GUIpics/File_menu.png" width="500" height="400"> | <img src="GUIpics/Help_menu.png" width="500" height="400"> |

If we click on the option Select Folder from the Folder/File Opts, a window will appear for selecting the desired folder with the PDB input files. Once specified, on the left frame we will see a list of those files. 

| *Selecting the input folder* | *List with input files* |
| :---: | :---: |
| <img src="GUIpics/selection_folder.png" width="500" height="400"> | <img src="GUIpics/list_files.png" width="500" height="400"> |

Then sepecifiying the desired options, we can click on the "Show options" button to see a summary of the different options. 

| *Specified options* | *Messagebox* |
| :---: | :---: |
| <img src="GUIpics/options.png" width="500" height="400"> | <img src="GUIpics/messagebox.png" width="400" height="200"> |

### Next steps

On the following weeks we intend to finish the GUI interface. For the moment when clicking on the "Run program!" button the program seems to run, but gives a weird error message from the Biopython module. 

```bash
/anaconda3/lib/python3.6/site-packages/Bio/PDB/Atom.py:125: RuntimeWarning: invalid value encountered in sqrt
  return numpy.sqrt(numpy.dot(diff, diff))
```
We're trying to fix the problem and have the application working as soon as possible!

## Requirements

In order to run this program with all its functionalities the user must have several packages dowloaded and working:

* [Python 3.6](https://www.python.org/downloads/)
* Modules: 
  * [Modeller v.9.19](https://salilab.org/modeller/download_installation.html)
  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)
  * [Biopython](http://biopython.org/wiki/Download)
  * argparse
  * os
  * sys 
  * numpy 
  

## Limitations

* Unable to model protein-RNA or protein-DNA interactions.
* Scalability: after certain recursion depth the program does not scalate too well. To improve the program, we could have used a dynammic programming approach, where out of several complexes only the most optimal is kept, and the next solutions are built onto that one.
* Stoichiometry: as powerful as it may be, the application relies too heavily on the stoichiometry input. On the one hand, this gives a greater control to the user over the output he or she expects but, on the other hand, it does not find complexes blindly, regardless of the stoichiometry.
* Without graphical interface for the moment. 
