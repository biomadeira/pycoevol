###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################
#TODO: 
# Interaction maps

import sys
from src import MAIN

def printUsage():
    "Prints the usage"
    __version__ = "beta"
    
    Usage = \
    """
    Pycoevol_%s (c) 2012, F. Madeira
  
    Pycoevol: A Python workflow to study protein-protein coevolution 
    and interaction.

    Pycoevol.py   -input1 -input2 -psiblast -alignment -coevolution
       
    input1        -seq1.fasta (-seqID1), -pdb1.pdb:A (-PDBID1:A)   
                  or -align1.fasta (where A is the chain designator)                
    input2        -seq2.fasta (-seqID2), -pdb2.pdb:B (-PDBID2:B)  
                  or -align2.fasta (where B is the chain designator) 
    psiblast      -internet, -local or -custom (NCBI's PSIBLAST and 
                  local database are optional) 
    alignment     -clustalw, -muscle, -mafft or -custom (MUSCLE and  
                  MAFFT are optional) 
    coevolution   -mi, -mie, -rcwmi,-cpvnmie, -cpvn, -clm, -vol
                  -omes, -pearson, -spearman, -mcbasc, -quartets,
                  -sca or -elsc
    help          -h or -help
                 
    Check the README.md for further details.
    """ %__version__
    print Usage
        
def pycoevolRun():
    "Routine which chooses the proper scripts given the input commands"
    main = MAIN.main(file1, file2,id1, id2, chain1, chain2, 
                     psiblast, alignment, coevolution)
    
    if psiblast == "custom" and alignment == "custom":
        main.coevolutionSripts()
    else:
        main.sequenceSripts()
        main.psiblastSripts()
        main.organismSripts()
        main.alignmentSripts()
        main.coevolutionSripts()
        main.infoScripts()
        
    print "Analysis complete!"
    return
   
def checkArguments():
    "Checks if the input commands are valid"
    try:
        input =str("./Data/" + file1)
        file = open(input,"r")
        file.close()
    except:
        #raise StandardError, "ERROR: File no.1 is not acessible"
        pass
    
    try:
        input =str("./Data/" + file2)
        file = open(input,"r")
        file.close() 
    except:
        #raise StandardError, "ERROR: File no.2 is not acessible"
        pass
    
    if psiblast != 'internet' and psiblast != 'local' and psiblast != 'custom':
        raise StandardError, "ERROR: PSI-Blast: Type 'internet', 'local'\
                            or 'custom'"
    
    if alignment != "clustalw" and alignment != "muscle" and \
    alignment != "mafft" and alignment != 'custom':
        raise StandardError, "ERROR: Alignment Tools: Type '-clustalw', \
        '-muscle', '-mafft' or 'custom'"
    
    if coevolution != 'mi' and coevolution != 'mie' and \
    coevolution != 'rcwmi' and coevolution != 'cpvnmie' and \
    coevolution != 'cpvn' and coevolution != 'clm' and \
    coevolution != 'vol' and coevolution != 'omes' and \
    coevolution != 'pearson' and coevolution != 'spearman' and \
    coevolution != 'mcbasc' and coevolution != 'quartets' and \
    coevolution != 'sca' and coevolution != 'elsc':
        raise StandardError, "ERROR: Coevolution Measure: Type '–mi', '–mie', \
        '–rcwmi', '–cpvnmie', '–cpvn', '–clm', '–vol', '-omes', '-pearson', \
        '-spearman', '-mcbasc', '-quartets', '-sca' or '-elsc'"

def checkDependencies():
    try: 
        import Bio
        del Bio
    except ImportError:
        raise ImportError, "ERROR: Unable to import Biopython"
    
    try: 
        import numpy
        del numpy
    except ImportError:
        raise ImportError, "ERROR: Unable to import Numpy"
    
    try: 
        import matplotlib
        del matplotlib
    except ImportError:
        raise ImportError, "ERROR: Unable to import Matplotlib"

def addtoPATH():
    sys.path.append("./src/tools")
    sys.path.append("./src/tools/blast+")
    sys.path.append("./src/tools/blast+/db")
    sys.path.append("./src/tools/clustalw")
    sys.path.append("./src/tools/mafft")
    sys.path.append("./src/tools/muscle")
    sys.path.append("./src/tools/phylip")
    
def main():
    "Pycoevol main program parses files and options as arguments"
    args = sys.argv
    
    global file1
    global chain1
    global id1
    global file2
    global chain2
    global id2
    global psiblast
    global alignment
    global coevolution
    
    if len(args) != 6:
        if args[0] == "-h" or args[0] == "-help":
            printUsage()
        else: 
            printUsage()
        return
    else:
        arg1 = args[1].split(":")
        if len(arg1) != 1: 
            file1 = arg1[0].lstrip("-")
            chain1 = arg1[1]
            id1 = file1.split(".")
            if len(id1) != 1:
                id1 = id1[0]
            else:
                id1 = id1[0]
                file1 = file1 + ".pdb"
        else:
            arg1 = args[1]
            file1 = arg1.lstrip("-")
            chain1 = ""
            id1 = file1.split(".")
            if len(id1) != 1:
                id1 = id1[0]
            else:
                id1 = id1[0] 
                file1 = file1 + ".fasta"
            
        arg2 = args[2].split(":")
        if len(arg2) != 1: 
            file2 = arg2[0].lstrip("-")
            chain2 = arg2[1]
            id2 = file2.split(".")
            if len(id2) != 1:
                id2 = id2[0]
            else:
                id2 = id2[0]
                file2 = file2 + ".pdb"
        else:
            arg2 = args[2]
            file2 = arg2.lstrip("-")
            chain2 = ""
            id2 = file2.split(".")
            if len(id2) != 1:
                id2 = id2[0]
            else:
                id2 = id2[0] 
                file2 = file2 + ".fasta"
        
        psiblast = args[3].lstrip("-")
        alignment = args[4].lstrip("-")
        coevolution = args[5].lstrip("-")
        
        checkArguments()
        checkDependencies()
        addtoPATH()
        pycoevolRun()
        return
    
if __name__ == "__main__":
    main()
    