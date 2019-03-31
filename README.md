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
- [Analysis of examples](#Analysis-of-examples)
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

Below is shown how to install and use this program as a stand-alone command line script (executing the script *MacroB.py*) or with the Graphical User Interface (*MB_GUI.py*).


### Software Requirements

These are the software and its versions required for the 4SMacroBuilder functionality and execution:

  * [Python 3.6](https://www.python.org/downloads/)
  * [Pymol](https://pymol.org/2/)


For the GUI the following ones are also necessary:

  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)



### Package tree

    Macrocomplex-Builder/
      README.md
      README.pdf
      setup.py
      MacroBuilder/
          __init__.py
          MacroB.py
          CustomPDB.py
          MB_GUI.py
          MB_launcher.py
      Examples/
          enterovirus/
          hemo/
          microtubul/
          nucleosome/
          phosphate/
          proteasome/
      templates/
          **all templates**
      doc/
          report.md
In order to be able to use all the scprits provided in 4SMacroBuilder the user has to install the package in the python site-packages.

```bash
   $ sudo python3 setup.py install
```
Be sure to have the dependencies previously stated.


### Installation via PIP

** To be discused **

### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, it can handle up to 95% of identity. 
* The same sequence appearing in different files with the same and/or different names. 
* Repeated chain interactions are not requiered as inputs (i.e. interaction A-A 10 times is treated as a single PDB). It solves infinite structures (i.e. Microtubul).


## Tutorial

In this section we make a brief explanation of how to use 4SMacroBuilder.

#### Command line Standalone application


```bash
    $ MBlauncher.py -h

    usage: MBlauncher.py [-h] -i, -input DIRECTORY [-o, -output OUTPUT]
                     [-c MAX_CHAINS] [-n NUM_MODELS] [-d] [-v] [-t TEMPLATE]

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
      -t TEMPLATE           To discriminate against different models, a template can be
                            given to calculate the RMSD
```

### GUI

Another way to use the program is using the the GUI. To do so run the following command:


```bash
$ MB_GUI.py
```
For a detailed explanation of how to use the GUI check the *report.pdf*

To get a better understanding of how to run the programme properly, we show different examples that represent different inputs that may be provided. The main aspects that may differ the inputs are: number of different chain interactions and number of atoms of the whole macrocomplex.


## Analysis of examples

Using the template (-t) optional argument. The program can compare between the model it's created and a given template. Using this we can analyse the quality of some reconstructed macrocomplexes.

### Proteosome

1pma PDB entry is a proteosome from *Thermoplasma acidophilum* (https://www.rcsb.org/structure/1PMA). A proteosome is a protein macrocomplex which degrade proteins by proteolysis of the protein peptide bonds. 1pma is a macrocomplex with two unique protein chains and a stoichiometry of hetero 28-mer-A14B14 with 4 interactions in one chain and 7 interactions in the other.
Throught our template analysis we do not see any differences with respect the number of chains between the one created and the template. Further analysis is done by structural visual comparison using an interactive visualization and molecular structural analysis programm such as UCSF CHIMERA. We can't oberve any differences in their structure, the superposition done is totally perfect. 

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/proteosome.png" alt="proteosome_superimposed_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 1</b></h4>
          <p><i>In blue we see the 1pma PDB protein and in brown the model created by 4SMacroBuilder</i></p>
        </div>
      </div>
    </div>
  </div>


### Enterovirus

3j23 PDB entry is the Enterovirus 71 empty capsids (https://www.rcsb.org/structure/3j23). EV71 os a single-stranded positive-sense RNA virus and a causative agent of hand, food, and mouth disease. 3j23 is a macrocomplex with three unique protein chains and a stoichiometry of hetero 180-mer-A60B60C60 with 4, 5 and 7 in each chain respectively.
Giving a set of protein-protein interactions, 4SMacroBuilder is able to construct the whole capsid macrocomplex. Comparing the structural composition of the model versus a template, we can observe any differences between both(image ***XX AND YY***).   

<div class="row">
    <div class="col-md-6">
      <div class="thumbnail">
        <img src="/images/3j23_enterovirus.png" alt="template_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 2</b></h4>
          <p><i>3j23 PDB Enterovirus 71 empty capsid</i></p>
        </div>
      </div>
    </div>
    <div class="col-md-6">
      <div class="thumbnail">
        <img src="/images/enterovirus.png" alt="enterovirus_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 3</b></h4>
          <p><i>Enterovirus model created by 4SMacroBuilder</i></p>
        </div>
      </div>
    </div>
</div>

### Nucleosome
3kuy PDB entry is the DNA stretching in the nucleosome core of *Escherichia coli* (https://www.rcsb.org/structure/3kuy).The DNA stretching in the nucleosome core can cause dramatic structural distortions, which may influence compaction and factor recognition in chromatin. It has a Stoichiometry of hetero 8-mer-A6B2.
4SMacroBuilder is able to create this protein - nucleic acid macrocomplex with 4 protein chains and 2 nucleotic acid chains with a total of 28 different pairwise interactions.

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/nucleosome.png" alt="model_nucleosome_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 4</b></h4>
          <p><i>Proteosome model created by 4SMacroBuilder</i></p>
        </div>
      </div>
    </div>
  </div>

## Strong Points

*1*. **Dynamic programming implementation**

The algorithm is based in a dynamic programming implementaton, in such a way that the final output is retrived in a very short time. This is due to the fact tha we use the interactions between chains as a previous knowlege to solve the final problem (*see documentation*). 

*2*. **Input managing**

  - The input names does not affect to the output (i.e. if all PDB files are named XY.pdb).
  - The input does not need all the interactions in different PDB files (i.e. case of virus capside or microtubul, with more than 150 chain-interaction in the case of the virus capside, and infinite interactions in the microtubul).
  - If the user gives a non-existing or wrong interaction the program ignores it and keeps going.

*3*. **Obtain different models**

Possibility to generate different models in a very short time. This allow the user to compare each model and decide which is/are the best models. Different models are scored in an output file.

*4*. **Launching the program with GUI**

4SMAcroBuilder can be launched from command line or with the **Graphical User Interface** (GUI). Besides, the GUI offer the advantatge to obtain a Pymol image of the final model, without the requirement of opening pymol.

*5*. **DNA & RNA interactions**
 
Possibility to model DNA/DNA, RNA/RNA, DNA/RNA, DNA/protein and RNA/protein interactions and retrieve a quick output (i.e. when modeing the ribosome).

*6*. **Modifiable number of chains in the final model**

Possibility to limit the number of chains when executing the program. Besides, if the user specifies that wants the macrocomplex with 7 chains but in fact the model has only 4 chains (i.e. Hemoglobin), it will not try to put more just because it was asked. This limited and reduces very much the program performance time.


## Limitations

*1*. **Increase of the computational cost with number of atoms in macrocomplex**

One of the main limitations dealing with the creation of a macrocomplex is the number of atoms and number of interactions it has. That's why we did a deeper anaylisis of these two factors using the microtubul folder. What is advantatgeous about this macrocomplex is that without any limitations it can go on forever without stopping, more or less like in a cell. But, limiting its parameters, it enables us to analyse our program.

We did a series of test normalizing by number of atoms and interactions. The microtuble has two different chains, with an average of 3347 atoms and 4 interactions by chain. 

As it can be seen in the following graph, the program follows an exponential curve. The more atoms/iterations it has to check, the more time it needs to run.

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/analysis.png" alt="time_analysis_image" style="width:600px;height:500px">
        <div class="caption">
          <h4><b>Figure 5</b></h4>
          <p><i>As we can see, the residues indicated by the arrow are *hydrophobic*, while the rest are hydrophilic.</i></p>
        </div>
      </div>
    </div>
  </div>

*2*. **Microtubul modeling**

Another factor that limit our program is that, due to some aspects of our approach, some "infinite" structures as the microtubule is not modeled as expected. This is possibly due to a random behavior implemented in the algorithm when adding subunits to the macrocomplex. 

*3*. **Different solutions**

Although the program can be asked to build more than one model from the same input, it is not able to deduce and build more than one output when there could be more than one possible solution. 

## Next Steps

*What could be the next future improvements?*

*1*. **Model Energy Minimization**

It could be implemented an option of energy optimization to a local energy minimum by molecular dynamics once the model/s has been finished. Protein structures often have errors of various magnitude: atoms overlapping, wrong side chain orientation (lack of water molecules when modeling). An energy minimization would look for the pathway that reduces the most in the overall energy of the system, obtaining a best approach of the final structure if possible. It could be reached with programs as **Amber**, that combines molecular mechanical force fields with biomolecular simularions. 
 
*2*. **Microtubul modeling**

It would be a good point to modify the algorithm approach which could improve the correct shape of the microtubul, as well as other non limit structures. We think that a way to do it could be to fisrt itearte the program by chain interactions as it does, but and, at a certain time force it to start again, but adding those interactions that had not been added yet.  

## FAQS

**Do I need PyMOL to launch and use the GUI?**

*Answer*:No, it is not necessary, but it is advisable to install [Pymol](https://pymol.org/2/) in order to see the whole macrocomplex/es once they have been done.