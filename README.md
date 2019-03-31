# MacrocomplexBuilder
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [MacrocomplexBuilder](#macrocomplexbuilder)
- [Software Requirements](#software-requirements)
- [Download and Installation](#download-and-installation)
  - [Package tree](#package-tree)
- [Input Files](#input-files)
- [Tutorial](#tutorial)
  - [Command line arguments](#command-line-arguments) 
  - [GUI](#gui)
- [Analysis of examples](#analysis-of-examples)
  - [Enterovirus](#enterovirus)
  - [Proteosome](#proteosome)
  - [Nucleosome](#nucleosome) 
- [Strong Points](#strong-points)
- [Computational Cost](#computational-cost)
- [Limitations](#limitations) 
- [Next steps](#next-steps)
- [FAQS](#FAQS)
  
<!-- /TOC -->



### MacrocomplexBuilder

MacrocomplexBuilder is a stand-alone python3 program developed by Pau Badia i Mompel, Altaïr C. Hernández and Natàlia Segura Alabart. It builds protein macrocomplexes taking a set of protein-protein, protein-RNA, protein-DNA, RNA - DNA, RNA - RNA, and/or DNA - DNA interactions. This software could serve to study quaternary structures that are difficult to study *in vivo*.

Below is shown how to install and use this program as a stand-alone command line script (executing the script *MacroB.py*) or with the Graphical User Interface (*MB_GUI.py*).


### Software Requirements

These are the software and its versions required for the MacrocomplexBuilder functionality and execution:

  * [Python 3.6](https://www.python.org/downloads/)
 

For the GUI the following ones are also necessary:

  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)
  * [Pymol](https://pymol.org/2/)


### Download and Installation

You can download our package using Git with the next command:
 
```bash
  $ git clone https://github.com/NataliaSegura/Macrocomplex-builder.git
  $ cd Macrocomplex-Builder
 ```
At this pont, the directory Macrocomplex-Builder should contain the files and directories described bellow:

#### Package tree

The package has the following structure:

    Macrocomplex-Builder/
      README.md
      README.pdf
      setup.py
      MB/
          __init__.py
          MacroB.py
          CustomPDB.py
          MB_GUI.py
          MB_launcher.py
      Examples/
          enterovirus/
          hemo/
          microtuble/
          nucleosome/
          phosphate/
          proteasome/
      Templates/
          3j23_enterovirus/
          1a3n_hemo/
          5syg_microtuble/
          3kuy_nucleosome/
          2f1d_phosphate/
          1pma_proteosome/
      doc/
          report.md

* README.md, README.pdf: the files containing the tutorial and information about our application.
* MB: a fold with the following scripts:
  - MBlauncher.py: the command-line script to launch the program.
  - MacroB.py: a module requiered by MBlauncher.py where are defined the classes of the program.
  - CustomPDB.py: a module required by MacroB.py where are defined the functions of the program.
  - MB_GUI.py: a module to launch the Graphical User Interface.
* Examples: a directory with several examples stored in sub-directories that serve as input to the program.
* Models: an empty folder where the created complexes will be saved.
* Templates: the raw PDB files from which we extracted the example pairwise interactions.
* setup.py script: to install the program in the python side-packages.
* doc: a folder with the report.md fie and report.pdf, as well as the documentation.

Check that all this information has been correctly downloaded and that there is the script called *setup.py*.

In order to be able to use all the scprits provided in MacrocomplexBuilder the user has to install the package in the python site-packages.

```bash
   $ sudo python3 setup.py install
```
Be sure to have the dependencies previously stated.


### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, it can handle 95% of identity. 
* The same sequence appearing in different files with the same and/or different names. 
* Repeated chain interactions are not requiered as inputs (i.e. interaction A-A 10 times is treated as a single PDB). It solves infinite structures (i.e. *Microtuble*).
* Pairwise interactions wrongly given to the program. The program threshold for considering two chains as interacting together is 3.5 Amstrongs. If the user gives interactions with bigger distance, they are not considered as such.
* A template PDB file containing the structure of the model to use it as a guideline.


All the models that only have protein - protein interactions can be dispalyed with both USCF CHIMERA and Pymol but when the model has nucleic acid interactions it must be opened with Pymol.


### Tutorial

In this section we make a brief explanation of how to use MacrocomplexBuilder.

#### Command line arguments

```bash
    $ MBlauncher.py -h

    usage: MBlauncher.py [-h] -i DIRECTORY [-o OUTPUT] [-c MAX_CHAINS]
                         [-n NUM_MODELS] [-d] [-v]
                         [-t TEMPLATE | -s STOICH_STRING]

    MacrocomplexBuilder is a python program designed to generate macrocomplex
    structures from simple pair inetractions

    optional arguments:
      -h, --help        show this help message and exit
      -i DIRECTORY      Input directory where the pair interaction pdbs are
                        located
      -o OUTPUT         Name of the output file, no extension is needed
      -c MAX_CHAINS     Maximum number of chains that the user wants in the model
      -n NUM_MODELS     Number of models that the program will compute
      -d, --dirty       Generates an output file for each added chain to track how
                        the program builds the complex
      -v, --verbose     Shows what the program is doing
      -t TEMPLATE       To discriminate against different models, a template can
                        be given to calculate the RMSD
      -s STOICH_STRING  The user can also give the desired stechiometry in this
                        format: A:6,B:11,C:2 ...
```

#### GUI

Another way to use the program is using the the GUI. To do so run the following command:

```bash
$ MB_GUI.py
```
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/whole_GUI.png" alt="whole_GUI" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>GUI structure</b></h4>
          <p>Here is the main structure of the GUI. First there is the Options widged where the user can select what parammeters to use for the modeling. Then there is the console panel where both the stderr and stdout will be shown, if the user wants more information there is the option verbose which will print each action that the program does. Next there is the Sequence widget, where the model sequences and their id's will be shown. At its bottom, there's the Structure composition panel where the model's chain composition is shown. Finally, a n image of the resulting model is shown at the bottom. (For Sequence, Structure and Image, all of the information shown is from the first model). </p>
        </div>
      </div>
    </div>
  </div>

To use the program, first a directory must be selected. One can do this by going to File>Select Directory:
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/selecting_directory.png" alt="select_dir" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Selecting a directory</b></h4>
          <p>To correctly select an input directory, the user must enter INSIDE the directory in the navigation window and press OK. </p>
        </div>
      </div>
    </div>
  </div>

Then, the user can change any of the options in the option panel, let's make just one model of a Proteosome to show it really works and select the option verbose to see what the program is doing. We press RUN to tun the program:
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/options_before_run.png" alt="options_b4_run" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Modify options before run</b></h4>
          <p>Here the model will make just one model, with the name proteasoma_1.cif and without checking structure composition. </p>
        </div>
      </div>
    </div>
  </div>

Finally we can see how the program has build this proteasoma really fast. 
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/final_result.png" alt="final_result" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Final result</b></h4>
          <p>The sequence, structure composition, structure and console have been updated</p>
        </div>
      </div>
    </div>
  </div>
### Analysis of examples

To get a better understanding of how to run the programme properly, we show different examples that represent different inputs that may be provided. The main aspects that may differ the inputs are: number of different chain interactions and number of atoms of the whole macrocomplex.

#### Enterovirus

The 3j23 PDB entry is the Enterovirus 71 empty capsids (https://www.rcsb.org/structure/3j23). EV71 is a single-stranded positive-sense RNA virus and a causative agent of hand, food, and mouth disease. 3j23 is a macrocomplex with three unique protein chains and a stoichiometry of hetero 180-mer-A60B60C60 with 4, 5 and 7 interaction sites in each chain respectively.
Giving a set of protein-protein interactions, MacrocomplexBuilder is able to construct the whole capsid macrocomplex. Comparing the structural composition of the model versus a template, we can observe any differences between both (*Figures 1*).   

To achieve this complex we can run:

```python
python3 MBlauncher.py -i enterovirus/ -o enterovirus -v
```

Where enterovirus/ is the Directory containing all input files, enterovirus is the file where the output will be saved in the current directory, and -v means that the standard error will be printed. 
This is a clear example of one of the strong points of our program: given 8 interaction pairs is able to produce the whole capsid with 180 chains in correct position.


<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/ent.png" alt="enterovirus_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 1</b></h4>
          <p><i>Comparison between the 3j23 PDB entry (left) and the model created by MacrocomplexBuilder(right)</i></p>
        </div>
      </div>
    </div>
  </div>
#### Proteosome

The 1pma PDB entry is a proteosome from *Thermoplasma acidophilum* (https://www.rcsb.org/structure/1PMA). A proteosome is a protein macrocomplex which degrade proteins by proteolysis of the peptide bonds. 1pma is a macrocomplex with two unique protein chains and a stoichiometry of hetero 28-mer-A14B14 with 4 interactions in one chain and 7 interactions in the other.

To achieve this complex we can run:

```python
python3 MBlauncher.py -i proteasome/ -o proteasoma -v -c 28
```

Where proteasome/ is the Directory containing all input files, proteasome is the file where the output will be saved in the current directory, -v means that the standard error will be printed, and 28 means that we are limiting the number of chains the model will have.

Even if in these case it is not necessary to limit the number of chains, we limited to 28 to show that the program will correctly construct the model. This will work with all the models this optional argument is given.

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/proteosome.png" alt="proteosome_superimposed_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 2</b></h4>
          <p><i>In blue we see the 1pma PDB protein and in brown the model created by MacrocomplexBuilder</i></p>
        </div>
      </div>
    </div>
  </div>


#### Nucleosome

As an optional argument, MacrocomplexBuilder can accept global stechiometry. If the user desires, it can be given and the program will create the model according to the given stechiometry.To achieve this complex we can run:

```python
python3 MBlauncher.py -i nucl/ -o nucl -v -s A6:B2
```
Where nucl/ is the Directory containing all input files, nucl is the file where the output will be saved in the current directory, -v means that the standard error will be printed, and A6:B2 means that the global stechiometry will be this one.

The 3kuy PDB entry is the DNA stretching in the nucleosome core of *Escherichia coli* (https://www.rcsb.org/structure/3kuy).The DNA stretching in the nucleosome core can cause dramatic structural distortions, which may influence compaction and factor recognition in chromatin. It has a Stoichiometry of hetero 8-mer-A6B2.
MacrocomplexBuilder is able to create this protein - nucleic acid macrocomplex with 4 protein chains and 2 nucleotic acid chains with a total of 28 different pairwise interactions.


<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/nucl.png" alt="model_nucleosome_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 3</b></h4>
          <p><i>Comparison between the nucleosome created limiting it stechiometry to A:6,B:2 (right) and the one created without stechiometry limitations (left)</i></p>
        </div>
      </div>
    </div>
</div>

### Strong Points

*1*. **Dynamic programming implementation**

The algorithm is based in a dynamic programming implementaton, in such a way that the final output is retrived in a very short time. This is due to the fact tha we use the interactions between chains as a previous knowlege to solve the final problem (*see documentation*). 

*2*. **Input managing**

  - The input names does not affect to the output (i.e. if all PDB files are named XY.pdb).
  - The input does not need all the interactions in different PDB files (i.e. case of virus capside or microtubul, with more than 150 chain-interaction in the case of the virus capside, and infinite interactions in the microtubul).
  - If the user gives a non-existing or wrong interaction the program ignores it and keeps going.

*3*. **Obtain different models**

Possibility to generate different models in a very short time. This allow the user to compare each model and decide which is/are the best models. Different models are scored in an output file.

*4*. **Launching the program with GUI**

MacrocomplexBuilder can be launched from command line or with the **Graphical User Interface** (GUI). Besides, the GUI offer the advantatge to obtain a Pymol image of the final model, without the requirement of opening pymol.

*5*. **DNA & RNA interactions**
 
Possibility to model DNA/DNA, RNA/RNA, DNA/RNA, DNA/protein and RNA/protein interactions and retrieve a quick output (i.e. when modeing the ribosome).

*6*. **Modifiable number of chains in the final model**

Possibility to limit the number of chains when executing the program. Besides, if the user specifies that wants the macrocomplex with 7 chains but in fact the model has only 4 chains (i.e. Hemoglobin), it will not try to put more just because it was asked. This limited and reduces very much the program performance time.

*7*. **Heteroatoms and water matter**

The active site of a protein often is composed by anions and cations. This information is described in the heteroatoms. MacrocomplexBuilder can use the heteroatom and water coordinates and information to construct the macrocomplex so we are not losing biological information. 

## Computational cost

One of the main limitations dealing with the creation of a macrocomplex is the number of atoms and number of interactions it has. That's why we did a deeper anaylisis of these two factors using the microtuble folder. What is advantatgeous about this macrocomplex is that without any limitations it can go on forever without stopping, more or less like in a cell. But, limiting its parameters, it enables us to analyse our program.

We did a series of tests normalizing by number of atoms and interactions. The microtuble has two different chains, with an average of 3347 atoms and 4 interactions by chain. 

Thanks to the illimited number of chains input that we can test when creating the microtuble we can asses its growing time curve. As it can be seen in the following graph, it might seem that until 200 chains the program followed a linear tendency but when a bigger number of calculations and steps was needed to create the microtuble, i.e. more chains, this behaviour is proven wrong. In fact, the program really follows an exponential curve. The more atoms/iterations it has to check, the more time it needs to run in an exopenential way.

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/analysis.png" alt="time_analysis_image" style="width:600px;height:500px">
        <div class="caption">
          <h4><b>Figure 4</b></h4>
          <p><i>Time performing program analysis. In the X-axis there is the number of chains the user can provide with the -t optional argument and in the Y-axis there is the time in seconds MacrocomplexBuilder needs.</i></p>
        </div>
      </div>
    </div>
  </div>


## Limitations

*1*. **Increase of the computational cost with number of atoms in macrocomplex**

As it can be seen in the *Figure 4*, the programs follows an exponential curve. The more atoms and interactions it has to check, the longer it takes to process.

*2*. **Microtuble modeling**

Another factor that limit our program is that, due to some aspects of our approach, some "infinite" structures are not modeled as expected, like the microtuble. This is possibly due to a random behavior implemented in the algorithm when adding subunits to the macrocomplex. 

*3*. **Different solutions**

Although the program can be asked to build more than one model from the same input, it is not able to deduce and build more than one output when there could be more than one possible solution. 

*4*. **The ATP problem and global stechiometry**

The problem with the ATP synthase macrocomplex is the number of interactions it has. The program can't handle all of them to create it and the wrong model is built. A way to modify the algorithmic approach is by givin stechiometry into the programm. That way, we limit the interactions and we force the macrocomplex into a specific shape. This can be achived using the optional argument -s (stechiometry). We give to the program the global macrocomplex stechiometry and it will build it using this parameters. A clear disatvantage of it is that even with the correct stechiometry it doesn't construct the right way because we are forcing the right number of chains but not whose interactions are in each one.


```bash
python3 MBlauncher.py -i ATP/ -o ATP_ -s A:1,B:1,C:1,D:1,E:2,F:1,G:3,H:1,I:8,J:1,K:1,L:1
```
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/atp.png" alt="atp_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 5</b></h4>
          <p><i>Comparison of 5arg PDB ATP Synthase model (right) and the model created with stechiometry limitations (left)</i></p>
        </div>
      </div>
    </div>
</div>

## Next Steps

*What could be the next future improvements?*

*1*. **Model Energy Minimization**

It could be implemented an option of energy optimization to a local energy minimum by molecular dynamics once the model/s has been finished. Protein structures often have errors of various magnitude: atoms overlapping, wrong side chain orientation (lack of water molecules when modeling). An energy minimization would look for the pathway that reduces the most in the overall energy of the system, obtaining a best approach of the final structure if possible. It could be reached with programs as **Amber**, that combines molecular mechanical force fields with biomolecular simularions. 
 
*2*. **Microtuble modeling**

It would be a good point to modify the algorithm approach which could improve the correct shape of the microtubule, as well as other non limit structures. We think that a way to do it could be to fisrt itearte the program by chain interactions as it does, but and, at a certain time force it to start again, but adding those interactions that had not been added yet.  

*3*. **ATP Synthase modeling**

The problem with these macrocomplex is the number of interactions it has and the program can't handle all of them to create it. A way to modify the algorithm approach to be able to construct correctly these macrocomplex is by givin stechiometry into the programm. THat way, we limit the interactions and we force the macrocomplex into a specific shape. One way to do it could be that given a template, the program calculates the stechiometry and use it to create the model.

## FAQS

**Do I need PyMOL to launch and use the GUI?**

*Answer*:No, it is not necessary, but it is advisable to install [Pymol](https://pymol.org/2/) in order to see the whole macrocomplex/es once they have been done.
