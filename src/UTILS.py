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
        
aa_hydrofobic = ['A', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W']
aa_hydrofile = ['C', 'N', 'Q', 'S', 'T', 'Y']
aa_basic = ['H', 'K', 'R']
aa_acid = ['D', 'E']

aa_symbols =  {'ALA':'A','CYS':'C','ASP':'D',
         'GLU':'E','PHE':'F','GLY':'G',
         'HIS':'H','LYS':'K','ILE':'I',
         'LEU':'L','MET':'M','ASN':'N',
         'PRO':'P','GLN':'Q','ARG':'R',
         'SER':'S','THR':'T','VAL':'V',
         'TYR':'Y','TRP':'W','XXX':'X'}

# Laguerre method with water. Esque et al, 2010
volume = {'N' : 125.2, 'P': 122.1, 'Q': 148.1,
          'A': 88.2, 'R': 188.8, 'S': 95.5,
          'C': 113.3,'T': 118.4, 'D': 113.4,
          'E': 134.8,'V': 134.5, 'F': 192.0,
          'W': 227.3,'G': 65.3,  'H': 159.2,
          'Y': 197.6,'I': 157.7, 'K': 164.2,
          'L': 158.7,'M': 164.9} 
