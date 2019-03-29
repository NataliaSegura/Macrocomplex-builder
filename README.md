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


The aim of this project is to build a protein macrocomplex (quaternary structure) with a single previous knowledge: a set of know protein pairwise interacions. Let's suppose that we know that there are the following protein - protein interactions: 

<div class="row">
    <div class="col-md-12" align="center">
      <div class="thumbnail">
        <img src="/approach1.png" alt="approach1_image" style="width:600px;height:500px">
        <div class="caption">
          <h4><b>Figure 1</b></h4>
          <p><i>Here A,B,C are chains/subunits of the macrocomplexe and may differ in sequence and structure</i></p>
        </div>
      </div>
    </div>
  </div>


In protein macrocomplexes there are several chains that interact with more than one chain, allowing the rest of interactions to be done. We could start taking one of these pair interactions as a template (i.e. A-B), and then superimpose the rest of the interactions by protein superposition. We have to assume that at least one chain of the template interact with another subunit (in this example A-C). This way we could superpose those identical chains (A-A) and move the new pair interaction to the template. Therefore, we would obtain a resulting structure of three chains (in this example). If we repeat this process until all the simillar chains are superposed, then we would obtain the final macrocomplex.

In order to carry this out we should know the order in which the program would have to superpose these pair of interactions, to avoid clashes between chains or even to prevent the program to superimpose the same chain more than one time. Also, we should know how many iterations the program would have to achive in order to make the final structure. 

This could be solved in different ways. For instance, starting with one pair of interactions as template it could checked all the possible interactions in each iteration an see which candidate would satisfy the problem statement(exhaustive search algorithm). Although this approach would be simple to implement, it would have a computational cost proportional to the number of candidate solutions, which would tend to grow in an exponential way.  

We approached this limitation focusing on the interacting residues of each subunit. If we imagine the whole macrocomplex structure as a lego puzzle, then we could realize that each chain has some residues that are interacting with at least another chain (hydrophobic residues) and the rest of the residues that are exposed to the solvent environment (hydrophilic residues). 

<div class="row">
    <div class="col-md-12" align="center">
      <div class="thumbnail">
        <img src="/approach2.png" alt="approach2_image" style="width:600px;height:500px">
        <div class="caption">
          <h4><b>Figure 2</b></h4>
          <p><i>As we can see, the residues indicated by the arrow are *hydrophobic*, while the rest are hydrophilic.</i></p>
        </div>
      </div>
    </div>
  </div>


If we base on that premise, we could store those interacting residues for each chain, as well as the corresponding chain those residues are interacting with. This is possible as we have this informations in the PDB files, so that would be the first task to do. In that way we could force the program to check in each iteration/superimposition whether those residues are interacting or not. At the same time we would have to consider as feasible complexes those that has no clashes when superimposed, which means that the backbone of the model that is been superimopsed is not interacting with the rest of the complex already joined (this means a threshold distance of 2 Aº).   

Then, our program would start to structurally superimpose structures with at least two identical subunits (those that share a pairwise sequence identity >= 95%), but this time for each model the program knows how many interacting sites are in each protein, and even with which specific chain has to interact on those sites. 

Let's explain this with an example: 

Imagine that we have the same pairwise interactions as before explained:

  * A-B
  * A-C
  * B-D

And know, when we read the PDB files, for each chain we store the interacting residues that take place by looking for those residues that are in a distance no larger that 5 Aº. Besides, we store in the same set which is the model of the corresponding interaction, so that we know which is the model to superimpose later. 

For instance:

               A -> (1,2,3) --> model A-B
                 -> (15,16) --> model A-C

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
