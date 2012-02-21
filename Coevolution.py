###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from os import system

def main():
    # example:
    # "python Pycoevol.py -3DX6.pdb:A -A2TP.pdb:A -internet -clustalw -sca")
    file1 = "3DX6.pdb"
    file2 = "A2TP.pdb"
    chain1 = "A"
    chain2 = "A"
    psiblast = "internet"
    alignment = "clustalw"
    coevolution = "sca"
    
    cline = system("python Pycoevol.py -%s:%s -%s:%s -%s -%s -%s") \
        %(file1, chain1, file2, chain2, psiblast, alignment, coevolution)
    print cline    
    
if __name__ == "__main__":
    main()
