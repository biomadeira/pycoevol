﻿PYCOEVOL
========
A Python workflow to study protein-protein coevolution and interaction
 
Pycoevol is an integrated system for studying inter-protein coevolution and interaction.
It automates the identification of contact points between protein partners, extending the 
general coevolution workflow consisting of: homologous sequence search; multiple sequence 
alignment computation; and coevolution analysis; with an improved selection of organisms 
and contact prediction. 

It generates friendly output results: matrix of scores; histograms;
heat-maps; PyMOL scripts and interaction maps. Additional information for common web-services
can be retrieved from SIFTS. 

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
[MUSCLE](http://www.drive5.com/muscle/),
[MAFFT](http://mafft.cbrc.jp/alignment/software/) and
[SIFTS lst files](http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)


##Usage
 
_python Pycoevol.py input1 input2 [options]_


 
	-h, --help		show this help message and exit
     
	-b PSIBLAST, --psiblast=PSIBLAST
 
					internet, local or custom
     
	-a ALIGNMENT, --alignment=ALIGNMENT
 
					clustalw, muscle, mafft or custom
     
	-c COEVOLUTION, --coevolution=COEVOLUTION
 
					mi, mie, rcwmi, cpvn, clm, vol, omes, pearson,spearman, mcbasc, quartets, sca or elsc
     
	-i IDS, --id=IDS
 
	-x CHAINS, --chain=CHAINS
 
	-p PARAMETERFILE, --parameters=PARAMETERFILE

For a detailed overview on how to install and use Pycoevol, please refer to the [User Guide] (https://sites.google.com/site/fmadeirawiki/home/pycoevol/Pycoevol_userguide.pdf) 


**Coevolution measures:**

* Mutual Information (mi) [Gloor et al, 2005]
* MI by pair Entropy (mie) [Martin et al, 2005]
* Row and Column Weighed MI (rcwmi) [Gouveia-Oliveira et al, 2007]
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

Source code available at
https://github.com/fmadeira/pycoevol


*Fábio Madeira and Ludwig Krippahl, 2012*