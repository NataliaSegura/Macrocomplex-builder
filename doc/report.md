# MacrocomplexBuilder 
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Introduction](#Introduction)
- [Background and Scientific explanation](#Background-and-Scientific-explanation)
- [Algorithm implementation](#Algorithm-implementation)
- [References](#References)
<!-- /TOC -->



### Introduction

Proteins are the most versatil macromolecules in living systems, as are in charge on multitude of specific chemical trasnformations, which provide the cell with usable energy and the molecules needed to form its structure and maintain the intracellular *homeostasis*. Proteins also recieve signals from outside the cell, starting intracellular signal transductions and regulating the gene expression in different stress situations. Sometimes, in order to perform their fucntion, the interaction of several protein subunits is made to form polypeptide complexes. Those structures are known as the quaternary structure of a protein, in which each subunit is atomically stabilized by hydrogen bonds or disulfide bonds, as well as other non bonding interactions (electrostatic and Van der Waals). Some examples are the hemoglobin, the ATP syntase, the RNA polymerase or the ribosome. 

Nevertheless, understanding how proteins interact with others in the assembly process is not an easy task. For this reason, different research groups have developed methods that predict how this interaction may occur. In the Protein Data Bank (PDB) are stored a large set of proteins with known structure, after a process of x-ray crystallography or Nuclear Magnetic Ressonance (NMR), allowing us to study the molecular space and possible allosteric interactions. 

It is important to stand out that there are some proteins that are difficult to crystallize due to its molecular conformation or dimensionallity, such as virus capsides. In other cases, it could be found some other type of interactions like protein/DNA, protein/RNA, DNA/DNA and/or RNA/RNA interactions, and it could be an interesting feature to study in transcriptional or translational processes (it would be the case of ribosomes, for instance). 

The main scope of this project is to reconstruct protein macrocomplexes from individual protein pairwise interactions using Bioinformatic resources. In order to do so, we have developed a stand-alone application that reads a set of protein - protein interactions in PDB file format and reconstruct different multi-subunit complexes. In this report we explain the approach we implement in order to make the program efficient and biologically reasonable. 


## Background and Scientific explanation

Macrocomplexes are built with an specific spatial order of subunits interacting with other subunits, and are usually stabilized by an hydrophobic core. This means that proteins interact between them with residues that don't like water (hydrophobic), and as a counterpart expose to the solvent those other residues that do like water (hydrophilic).

If we start from scratch to reconstruct a complex of subunits we need to know at least the number of chains that will build that macrostructure. Let's suppose we have the following protein - protein interactions: 

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/approach.png" alt="approach_explanation_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 1</b></h4>
          <p><i>Here A,B,C are chains/subunits of the macrocomplex and may differ in sequence and structures. </i></p>
        </div>
      </div>
    </div>
  </div>

  
In protein macrocomplexes there are several chains that interact with more than one chain, allowing the rest of interactions to be done. We could start by taking one of these pair interactions as a template (i.e. A-B), and then superimpose the rest of the interactions by protein superposition. We have to assume that at least one chain of the template interact with another subunit (in this example A-C). This way we could superpose those identical chains (A-A) and move the new pair interaction to the template. Therefore, we would obtain a resulting structure of three chains (*Figure 1*). If we repeat this process until all the simillar chains are superposed, then we would obtain the final macrocomplex.

In order to carry this out we should know the order in which the program would have to superpose these pair of interactions. The order is needed to avoid clashes between chains or even to prevent the program to superimpose the same chain more than onece. Also, we should know how many iterations the program would have to achive in order to make the final structure. 

This could be solved in different ways. For instance, starting with one pair of interactions as template it could check all the possible interactions in each iteration an see which candidate would satisfy the problem statement. This approach is called exhaustive search algorithm. Although this approach would be the simplest and easiest to implement, it would have a computational cost proportional to the number of candidate solutions, which would tend to grow in exponentially.  


### Algorithm implementation


We approached the exhaustive search algorithm by focusing on the interacting residues of each subunit. If we imagine the whole macrocomplex structure as a lego puzzle, then we would realize that each chain has some residues that are interacting with at least another chain (*hydrophobic residues*) and the rest of the residues that are exposed to the solvent environment (*hydrophilic residues*). 

Basing on that premise, we stored the interacting residues for each chain, as well as the corresponding chain those residues are interacting with.  These takes place by looking for those residues that are in no larger that 5 Amstrongs distance. This is possible because we have this information in the PDB files, so that would be the first task to do. In that way we could force the program to check in each iteration/superimposition whether those residues are interacting or not. At the same time, we would have to consider as feasible complexes those that has no clashes when superimposed, which means that the backbone of the model that is been superimopsed is not interacting with the rest of the complex already joined - a threshold distance of 2 Amstrongs.   

Then, our program would start to structurally superimpose structures with at least two identical subunits (those that share a pairwise sequence identity bigger or equal to 95%). For each model the program knows how many interacting sites are in each protein, and even with which specific chain has to interact on those sites. 

And know, when we read the PDB files, for each chain we load into memory the interacting residues. Besides, we store in the same set which is the model of the corresponding interaction, so that we know which is the model to superimpose later. At the end we will have a dictionary with each chain as key and all the information as values: interacting residues, interacting chain pair names, interacting residues for the other chain. Once we have this we update the interaction dictionary of each chain inside the dictionary to have the information in it. One point to remark is that when the updating is done the same interaction is not inserted twice or more.
The superimposition will start with the chain with more interactions, to avoid starting for a wrongly given interaction that just interact with itself and once the initial model is chosen the rest will be randomly added. When the superposing occurs we use an internal checklist to know for each chain which interactions are done and which ones are to be checked and superimposed. But if a chain appears more than once in the complex, each of them will have its checklist to complete. That way we ensure that all interactions are done.

### Examples

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

### References

Since our project needed a way to determine if two residues where interacting or to check if a clash had been produced, we used the following values from the following references:

* In each interaction pair we store the atoms from each chain that interacts with the other. Since many PDB files lack hydorgen bonds, the presence of an hydrogen bond can be inferred when an atom from one chain is within 3.5 Amstrongs of an atom of the other chain.  
Martz, Eric; Help, Index & Glossary for Protein Explorer, http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/igloss.htm
Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.

* Before two chains are superimposed, in order to determinate if they will be or not, we will check if they have clashes. The threshold that we use is by measuring the atomic distance. If the distance is 2 Amstrongs or less and there is more than a 3% of the alpha carbons atoms in proteins and the carbon one in nucleic acid atoms has clashes we will not superpose them.  
Values of Vanderwalls radius taken from: http://ww2.chemistry.gatech.edu/~lw26/structure/molecular_interactions/mol_int.html
Batsanov S.S.; Van der Waals Raddi of Elements, Inorganic Materials, 2001.
