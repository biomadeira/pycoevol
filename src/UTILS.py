# Encoding utf-8
# Created by F. Madeira, 2012
# This code is part of Pycoevol distribution.
# This work is public domain. 
"""
Utilities used in some routines.
"""

aa = ['A','C','D','E','F','G','H', 'K','I','L','M','N','P','Q','R','S','T','V','Y','W']

aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'LYS', 'ILE', 'LEU',
           'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TYR', 'TRP']

aa_symbols =  {'ALA':'A','CYS':'C','ASP':'D',
         'GLU':'E','PHE':'F','GLY':'G',
         'HIS':'H','LYS':'K','ILE':'I',
         'LEU':'L','MET':'M','ASN':'N',
         'PRO':'P','GLN':'Q','ARG':'R',
         'SER':'S','THR':'T','VAL':'V',
         'TYR':'Y','TRP':'W','XXX':'X'}      
  
# amino acid properties      
aa_hydrofobic = ['A','F','G','I','L','M','P','V','W']
aa_hydrofile = ['C','N','Q','S','T','Y']
aa_basic = ['H','K','R']
aa_acid = ['D','E']
aa_polar = ['S','T','Q','C','E','Y','D','K','H','R','N']
aa_non_polar = ['A','V','L','I','G','W','F','P','M']
aa_charged = ['H','R','K','D','E']

# Hydrophobicity scale:
# Kyte J and Doolittle RF: A simple method for displaying the 
# hydropathic character of a protein. J Mol Biol 157:105, 1982.
kyte_doolittle = {'A':1.8,'C':2.5,'D':-3.5,
                  'E':-3.5,'F':2.8,'G':-0.4,
                  'H':-3.2,'K':-3.9,'I':4.5,
                  'L':3.8,'M':1.9,'N':-3.5,
                  'P':-1.6,'Q':-3.5,'R':-4.5,
                  'S':-0.8,'T':-0.7,'V':4.2,
                  'Y':-1.3,'W':-0.9}

# Hoop TP and Woods KR: Prediction of protein antigenic determinants 
# from amino acid sequences. Proc Natl Acad Sci USA 78:3824, 1981. 
hopp_woods = {'A':-0.5,'C':-1.0,'D':3.0,
              'E':3.0,'F':-2.5,'G':0.0,
              'H':-0.5,'K':3.0,'I':-1.8,
              'L':-1.8,'M':-1.3,'N':0.2,
              'P':0.0,'Q':0.2,'R':3.0,
              'S':0.3,'T':-0.4,'V':-1.5,
              'Y':-2.3,'W':-3.4}  

# D. Eisenberg; R. M. Weiss & T. C. Terwilliger:
# The hydrophobic moment detects periodicity in protein hydrophobicity.
# Proc Natl Acad Sci U S A, 81, 140-144
eisenberg = {'A':0.62,'C':0.29,'D':-0.9,
            'E':-0.74,'F':1.19,'G':0.48,
            'H':-0.4,'K':1.38,'I':-1.5,
            'L':1.06,'M':0.64,'N':-0.78,
            'P':0.12,'Q':-0.85,'R':-2.53,
            'S':-0.18,'T':-0.05,'V':1.08,
            'Y':0.81,'W':0.26}  

# D. M. Engelman; T. A. Steitz & A. Goldman:
# Identifying nonpolar transbilayer helices in amino acid sequences of 
# membrane proteins. Annu Rev Biophys Biophys Chem, 15, 321-353
engelman = {'A':1.6,'C':2.0,'D':-9.2,
            'E':-82,'F':3.7,'G':1.0,
            'H':-3.0,'K':3.1,'I':-8.8,
            'L':2.8,'M':3.4,'N':-4.8,
            'P':-0.2,'Q':-4.1,'R':-12.3,
            'S':0.6,'T':1.2,'V':2.6,
            'Y':1.9,'W':-0.7}

# J. L. Cornette; K. B. Cease; H. Margalit; J. L. Spouge; J. A. Berzofsky & C. DeLisi:
# Hydrophobicity scales and computational techniques for detecting amphipathic 
# structures in proteins. J Mol Biol, 195, 659-685
cornette = {'A':0.2,'C':4.1,'D':-3.1,
            'E':-1.8,'F':4.4,'G':0.0,
            'H':0.5,'K':4.8,'I':-3.1,
            'L':5.7,'M':4.2,'N':-0.5,
            'P':-2.2,'Q':-2.8,'R':1.4,
            'S':-0.5,'T':-1.9,'V':4.7,
            'Y':1.0,'W':3.2}


# Amino acid's volume:
# Laguerre method with water. Esque et al, 2010
volume = {'N' : 125.2, 'P': 122.1, 'Q': 148.1,
          'A': 88.2, 'R': 188.8, 'S': 95.5,
          'C': 113.3,'T': 118.4, 'D': 113.4,
          'E': 134.8,'V': 134.5, 'F': 192.0,
          'W': 227.3,'G': 65.3,  'H': 159.2,
          'Y': 197.6,'I': 157.7, 'K': 164.2,
          'L': 158.7,'M': 164.9} 

# atomic radius
radii = {'H': 1.20, 'N': 1.55, 'NA': 2.27,
         'CU': 1.40, 'CL': 1.75, 'C': 1.70,
         'O': 1.52, 'I': 1.98, 'P': 1.80,
         'B': 1.85, 'BR': 1.85, 'S': 1.80,
         'SE': 1.90, 'F': 1.47, 'FE': 1.80,
         'K':  2.75, 'MN': 1.73, 'MG': 1.73,
         'ZN': 1.39, 'HG': 1.8, 'XE': 1.8,
         'AU': 1.8, 'LI': 1.8, '.': 1.8}
