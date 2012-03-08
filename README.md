PYCOEVOL
========

A Python workflow to study protein-protein coevolution and interaction

##Disclaimer 

This software is provided "as is", with no explicit or implied warranties. 
Use this software at your own risk.

##Copyright

This software is public domain, and everyone has the right to copy, 
distribute, reuse, modify, improve and debug it.

If you want to cite this piece of software/workflow use the following:

Fábio Madeira and Ludwig Krippahl. 2012. PYCOEVOL: A Python workflow to study 
protein-protein coevolution. Proceedings of the International conference on 
Bioinformatics Models, Methods and Algorithms - BIOINFORMATICS 2012, pp.143-9. 

##Dependencies

[Python 2.7.2](http://python.org/),
[Biopython 1.58](http://biopython.org/),
[Numpy 1.6.1](http://numpy.scipy.org/),
[Matplotlib 1.1.0](http://matplotlib.sourceforge.net/) and
[ClustalW](http://www.clustal.org/)

**Optional:**
[NCBI Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download),
[NCBI's "refseq_protein" database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/),
[MUSCLE](http://www.drive5.com/muscle/) and
[MAFFT](http://mafft.cbrc.jp/alignment/software/)


##Usage
 
_python Pycoevol.py  -input1 -input2 -psiblast -alignment -coevolution_
    
    input1        -seq1.fasta (-seqID1), -pdb1.pdb:A (-PDBID1:A)   
                  or -align1.fasta (where A is the chain designator)                
    input2        -seq2.fasta (-seqID2), -pdb2.pdb:B (-PDBID2:B)  
                  or -align2.fasta (where B is the chain designator) 
    psiblast      -internet, -local or -none (NCBI's PSIBLAST and 
                  local database are optional) 
    alignment     -clustalw, -muscle, -mafft or -none (MUSCLE and  
                  MAFFT are optional) 
    coevolution   -mi, -mie, -rcwmi,-cpvnmie, -cpvn, -clm, -vol
                  -omes, -pearson, -spearman, -mcbasc, -quartets,
                  -sca or -elsc
    help          -h or -help

**Examples:**

_python Pycoevol.py -3DX6.pdb:A -3DX6.pdb:B -internet -clustalw -sca_

_python Pycoevol.py -3DX6:A -3DX6:B -internet -clustalw -sca_

_python Pycoevol.py -seq1.fasta -seq2.fasta -local -mafft -mi_

_python Pycoevol.py -P30481 -P61769 -local -mafft -mi_

Pycoevol.py is the main execution file which parses files and options
as arguments. Parameters.py must be edited in order to specify some 
parameters.

The processed files are written in the ./Pycoevol/Data folder and 
output results are written in the ./Pycoevol/Results folder.
Make sure you insert input files in the ./Pycoevol/Data folder.

**Coevolution measures:**

* Mutual Information (mi) [Gloor et al, 2005]
* MI by pair Entropy (mie) [Martin et al, 2005]
* Row and Column Weighed MI (rcwmi) [Gouveia-Oliveira et al, 2007]
* Contact Preferences, Volume Normalized MIE (cpvnmie) [F. Madeira, 2012 - unpublished]
* Contact Preferences, Volume Normalized (cpvn) [Glaser et al, 2001]
* Contact PDB-derived Likelihood Matrix (clm) [Singer et al, 2002]
* Residue-residue Volume Normalized (vol) [based on Esque et al, 2010]
* Observed Minus Expected Squared  (omes) [Kass and Horovitz, 2002]
* Pearson’s correlation (pearson) [Göbel et al, 1994]
* Spearman’s rank correlation (spearman) [Pazos et al, 1997]
* McLachlan Based Substitution Correlation (mcbasc) [Fodor and Aldrich, 2004]
* Quartets (quartets) [Galitsky, 2002]
* Statistical Coupling Analysis (sca) [Lockless and Ranganathan, 1999]
* Explicit Likelihood of Subset Covariation (elsc) [Dekker et al, 2004]

**Pairwise distance measures:**

* ClustalW distance[Chenna et al, 2003]
* p-distance [Jukes and Cantor, 1969]
* Jukes-Cantor [Jukes and Cantor, 1969]
* Kimura distance [Kimura, 1983]
* Pairwise score using Dayhoff or Henikoff matrices [Dayhoff et al, 1978; 
Henikoff and Henikoff, 1992]

##Availability

Packages available at
http://code.google.com/p/pycoevol/

Source code available at
https://github.com/fmadeira/pycoevol


*F. Madeira, 2012*

