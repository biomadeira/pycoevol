PYCOEVOL
========

Madeira, F. and Krippahl, L. PYCOEVOL: A Python workflow to study 
protein-protein coevolution. BIOINFORMATICS 2012 (accepted)

##Disclaimer 

This software is provided "as is", with no explicit or implied 
warranties. Use this software at your own risk.

##Copyright

This software is public domain, and everyone has the right to copy, 
distribute, reuse, modify, improve and debug it.

If you want to cite this piece of software/workflow use the following:
Madeira, F., Krippahl, L. (2011). PYCOEVOL: A Python workflow to study 
protein-protein coevolution. International conference on Bioinformatics 
Models, Methods and Algorithms - BIOINFORMATICS 2012, Vila Moura, 
Portugal, February 1-4, 2012 

##Dependencies

[Python 2.7.2](http://python.org/),
[Biopython 1.58](http://biopython.org/),
[Numpy 1.6.1](http://numpy.scipy.org/),
[Matplotlib 1.1.0](http://matplotlib.sourceforge.net/) and
[ClustalW](http://www.clustal.org/)

Optional:
[NCBI Blast+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
[NCBI "refseq_protein" database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/),
[Muscle](http://www.drive5.com/muscle/) and
[Mafft](http://mafft.cbrc.jp/alignment/software/)


##Usage
 
_python Pycoevol.py  -input1 -input2 -psiblast -alignment -coevolution_
       
    input1        -sequence1.fasta or -pdb1.pdb:A, where A is the 
                  chain designator                  
    input2        -sequence2.fasta or -pdb2.pdb:B, where B is the 
                  chain designator
    psiblast      -internet or -local (needs NCBI's PSIBLAST and 
                  local database)  
    alignment     -clustalw, -muscle or -mafft (MUSCLE and MAFFT 
                  are optional)
    coevolution   -mi, -mie, -rcwmi,-cpvnmie, -cpvn, -clm or -vol
    help          -h or -help

**Example:** 
_python Pycoevol.py -3DX6.pdb:A -A2TP.pdb:A -internet -clustalw -mi_

Pycoevol.py is the main execution file which parses files and options
as arguments. Parameters.py must be edited in order to specify some 
parameters.

The processed files are written in the ./Pycoevol/Data folder and 
output results are written in the ./Pycoevol/Results folder.
Make sure you insert input files in the ./Pycoevol/Data folder.

##Availability

Packages available at
http://code.google.com/p/pycoevol/

Source code available at
https://github.com/fmadeira/pycoevol


*F. Madeira, 2012*

