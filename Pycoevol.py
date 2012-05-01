###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################
#TODO: 
# Interaction maps

import os
import sys
from src import MAIN
from Parameters import LoadParameters as LP
from optparse import OptionParser
from Bio.Align.Applications import ClustalwCommandline

def printUsage():
    """Prints the usage - DEPRECATED"""
    __version__ = "beta"
    
    Usage = \
    """
    Pycoevol_%s (c) 2012, F. Madeira
  
    Pycoevol: A Python workflow to study protein-protein coevolution 
    and interaction.

    Pycoevol.py   input1 input2 
       
    input1        seq1.fasta (-seqID1), pdb1.pdb:A (-PDBID1:A)   
                  or align1.fasta (where A is the chain designator)                
    input2        seq2.fasta (seqID2), -pdb2.pdb:B (-PDBID2:B)  
                  or -align2.fasta (where B is the chain designator) 
    -p --psiblast
                  internet, local or custom (NCBI's PSIBLAST and 
                  local database are optional) 
    -a --alignment
                  clustalw, muscle, mafft or custom (MUSCLE and  
                  MAFFT are optional) 
    -c --coevolution
                  mi, mie, rcwmi, cpvnmie, cpvn, clm, vol
                  omes, pearson, spearman, mcbasc, quartets,
                  sca or elsc
    -x --chain
                  chain identifier (in same order as input file). Default A
    -i --id
                  identifier for each protein, in same order as input files.
    -h --help
                 
    Check the README.md for further details.
    """ % __version__
    print Usage
        
def pycoevolRun():
    "Routine which chooses the proper scripts given the input commands"
    main = MAIN.main(file1, file2, id1, id2, chain1, chain2, parameterfile,
                     psiblast, alignment, coevolution, dirname)
    
    if psiblast == "custom" and alignment == "custom":
        print 'Coevolution scripts...'
        sys.stdout.flush()
        main.coevolutionSripts()
        print '... OK'
    else:
        print 'Sequence scripts...'
        sys.stdout.flush()
        main.sequenceSripts()
        print '... OK'
        
        print 'BLAST scripts...'
        sys.stdout.flush()
        main.psiblastSripts()
        print '... OK'
        
        print 'Organism scripts...'
        sys.stdout.flush()
        main.organismSripts()
        print '... OK'
        
        print 'Alignment scripts...'
        sys.stdout.flush()
        main.alignmentSripts()
        print '... OK'
        
        print 'Coevolution scripts...'
        sys.stdout.flush()
        main.coevolutionSripts()
        print '... OK'
        
        print 'Info scripts...'
        sys.stdout.flush()
        main.infoScripts(SIFTS)
        print '... OK'
    return
   
def checkArguments():
    "Checks if the input commands are valid"
    try:
        input = str("./Data/" + file1)
        file = open(input, "r")
        file.close()
    except:
        #raise StandardError, "ERROR: File no.1 is not acessible"
        pass
    
    try:
        input = str("./Data/" + file2)
        file = open(input, "r")
        file.close() 
    except:
        #raise StandardError, "ERROR: File no.2 is not acessible"
        pass
    
    if len(chain1) <= 2 and len(chain2) <= 2:
        pass
    else:
        raise StandardError, "ERROR: Chains' length must be = 1"
    
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
    "Checks the import of mandatory python modules and clustalw"
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
    
    try:
        input = "./src/tools/clustalw/test/test.fasta"
        clustalw = ClustalwCommandline(infile=input) 
        clustalw()
        os.remove("./src/tools/clustalw/test/test.aln")
        os.remove("./src/tools/clustalw/test/test.dnd")
    except:
        raise StandardError, "ERROR: Unable to run ClustalW"
        

def checkSIFTS():
    "Checks the availability of SIFTS files"
    global SIFTS
    try:
        input = str("./SIFTS/pdb_chain_scop_uniprot.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_chain_cath_uniprot.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_chain_enzyme.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_chain_interpro.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_chain_pfam.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_chain_taxonomy.lst")
        file = open(input, "r")
        file.close()
        input = str("./SIFTS/pdb_pubmed.lst")
        file = open(input, "r")
        file.close()
        SIFTS = True
        print "SIFTS... OK"
    except:
        SIFTS = False        
        print "SIFTS... NOT OK"
        
def addtoPATH():
    sys.path.append("./src/tools")
    sys.path.append("./src/tools/blast+")
    sys.path.append("./src/tools/blast+/db")
    sys.path.append("./src/tools/clustalw")
    sys.path.append("./src/tools/mafft")
    sys.path.append("./src/tools/muscle")
    
def ParseArguments():
    global file1
    global id1
    global chain1
    global file2
    global id2
    global chain2
    global parameterfile
    global psiblast
    global alignment
    global coevolution
    global dirname

    # defaults
    pathcwd = os.getcwd()
    dirname = os.getcwd() + "/Results/"
    parameterfile = ''
    file1 = ''
    file2 = '' 
    chain1 = ''
    chain2 = ''

    parser = OptionParser(usage='Pycoevol.py input1 input2 [options]')
    parser.add_option('-b', '--psiblast', type='string',
                      dest='psiblast', default='internet',
                      help='internet, local or custom')
    parser.add_option('-a', '--alignment', type='string',
                      dest='alignment', default='clustalw',
                      help='clustalw, muscle, mafft or custom')
    parser.add_option('-c', '--coevolution', type='string',
                      dest='coevolution', default='mi',
                      help='mi, mie, rcwmi, cpvn, clm, vol, omes, pearson, spearman, mcbasc, quartets, sca or elsc')
    parser.add_option('-i', '--id', action='append', type='string',
                      dest='ids', default=[])
    parser.add_option('-x', '--chain', action='append', type='string',
                      dest='chains', default=[])
    parser.add_option('-p', '--parameters',
                      dest='parameterfile', default=parameterfile)
      
    (options, args) = parser.parse_args()
    if len(args) == 0 and len(options.ids) == 0:
        parser.print_help()
        sys.exit()
        
    if len(args) == 2:
        input1 = args[0]
        input2 = args[1]
        dirname = os.path.dirname(input1) + "/"
        file1 = os.path.basename(input1)
        file2 = os.path.basename(input2)
        id1 = file1.split(".")[0]
        id2 = file2.split(".")[0]
    if len(options.chains) == 2:
        chain1 = options.chains[0]
        chain2 = options.chains[1]
    if len(options.ids) == 2:
        id1 = options.ids[0]
        id2 = options.ids[1]
        if chain1 == '' and chain2 == '':
            file1 = id1 + ".fasta"
            file2 = id2 + ".fasta"
        else:
            file1 = id1 + ".pdb"
            file2 = id2 + ".pdb"
    if options.parameterfile != '':
        parameterfile = options.parameterfile.strip('"')
        LP(parameterfile, "test")
    else:
        parameterfile = pathcwd + "/Params.config"
        parameterfile = parameterfile.strip('"')
        LP(parameterfile, "test")
    psiblast = options.psiblast
    alignment = options.alignment
    coevolution = options.coevolution

def main():
    ParseArguments()        
    checkArguments()
    print 'Arguments... OK'
    checkDependencies()
    print 'Dependencies... OK'
    checkSIFTS()
    addtoPATH()        
    pycoevolRun()
    print 'Analysis Complete !!'
    return
    
if __name__ == "__main__":
    main()
    

