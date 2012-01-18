###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

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
                 
    Check the README.txt for further details.
    """ %__version__
    print Usage
        
def pycoevolRun():
    "Routine which chooses the proper scripts given the input commands"
    main = MAIN.main(file1,file2,id1, id2, chain1, chain2, 
                     psiblast, alignment, coevolution)
    main.sequenceSripts()
    main.psiblastSripts()
    main.organismSripts()
    main.alignmentSripts()
    main.coevolutionSripts()
    main.infoScripts()
    print "Analysis complet!!"
    
        
def checkArguments():
    "Checks if the input commands are valid"
    try:
        input =str("./Data/" + file1)
        file = open(input,"r")
        file.close()
    except:    
        raise StandardError, "ERROR: File no.1 is not acessible"
    
    try:
        input =str("./Data/" + file2)
        file = open(input,"r")
        file.close() 
    except:    
        raise StandardError, "ERROR: File no.2 is not acessible"
    
    if psiblast != 'internet' and psiblast != 'local':
        raise StandardError, "ERROR: PSI-Blast: Type 'internet' or 'local'"
    
    if alignment != "clustalw" and alignment != "muscle" and \
    alignment != "mafft":
        raise StandardError, "ERROR: Alignment Tools: Type '-clustalw', \
        '-muscle' or '-mafft'"
    
    if coevolution != 'mi' and coevolution != 'mie' and \
    coevolution != 'rcwmi' and coevolution != 'cpvnmie' and \
    coevolution != 'cpvn' and coevolution != 'clm' and coevolution != 'vol':
        raise StandardError, "ERROR: Coevultion Measure: Type '–mi', '–mie', \
        '–rcwmi', '–cpvnmie', '–cpvn', '–clm' or '–vol'"

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
            id1 = id1[0]
        else:
            arg1 = args[1]
            file1 = arg1.lstrip("-")
            chain1 = ""
            id1 = file1.split(".")
            id1 = id1[0]
            
        arg2 = args[2].split(":")
        if len(arg2) != 1:    
            file2 = arg2[0].lstrip("-")
            chain2 = arg2[1]
            id2 = file2.split(".")
            id2 = id2[0]
        else: 
            arg2 = args[2]
            file2 = arg2.lstrip("-")
            chain2 = ""
            id2 = file2.split(".")
            id2 = id2[0]
            
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
    