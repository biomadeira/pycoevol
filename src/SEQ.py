###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

import src.SASA as sasa
from src.UTILS import aa_list, aa_symbols
from Parameters import LoadParameters as LP
import time
from os import remove
from shutil import copyfile
from urllib import urlopen
from Bio import SeqIO, Entrez
from Bio.Alphabet import IUPAC
from Bio.PDB.PDBParser import PDBParser
Entrez.email = "entrez@mail.com"

class sequence:
    """
    Main code for handling sequences and structures.
    """
    def __init__(self, file1, file2, id1, id2, chain1, chain2, parameterfile, 
                 dirname):
        self.file1 = file1
        self.file2 = file2
        self.chain1 = chain1
        self.chain2 = chain2
        self.id1 = id1
        self.id2 = id2
        self.parameterfile = parameterfile
        self.dirname = dirname
        
    def __call__(self, file1, file2, id1, id2, chain1, chain2, parameterfile,
                 dirname):
        self.file1 = file1
        self.file2 = file2
        self.chain1 = chain1
        self.chain2 = chain2
        self.id1 = id1
        self.id2 = id2
        self.parameterfile = parameterfile
        self.dirname = dirname
    
    def validFASTA(self, file, id):
        "Checks if the input file is a valid FASTA file"
        
        try:
            input = str(self.dirname + file)
            SeqIO.read(input, "fasta", IUPAC.protein)
        except:
            try:
                "Fetches a sequence according to GI identifier or UniProt ID"
                fetch = Entrez.efetch(db="protein", id=id, rettype="fasta")
                              
                output = str(self.dirname + file)
                out = open(output, "w")
                out.write(fetch.read())
                out.close()
                read = SeqIO.parse(output, "fasta", IUPAC.protein)
                for record in read:
                    sequence = str(record.seq)
                out = open(output, "w")
                print >> out, ">Query_id" + "\n" + sequence + "\n"
                out.close()
                
            except:                           
                raise StandardError, "%s - Invalid sequence identifier or sequence file" % (id)                          
    
    def queryFASTA(self, file, id):
        "Changes FASTA original header to 'Query_id'"
        
        input = str(self.dirname + file)
        input_sequence = SeqIO.parse(input, "fasta", IUPAC.protein)
        for record in input_sequence:
            sequence = str(record.seq)
            break
            
        output = str(self.dirname + id + ".fa")
        out = open(output, "w")
        print >> out, ">Query_id" + "\n" + sequence + "\n"
        out.close()
        input_sequence.close()
        
        remove(input)
        
    def validPDB(self, file, id, chain):
        "Checks if the input file is a valid PDB file"
    
        try:
            input = str(self.dirname + file)
            PDBParser().get_structure(id, input)
            try:
                test_structure = PDBParser().get_structure(id, input)
                test_model = test_structure[0]
                test_model[chain]
            except: 
                raise StandardError, "%s - Invalid chain" % (chain)
        except: 
            try:
                "Fetches a PDB file from the RCSB Protein Databank"
                url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
                read = urlopen(url).read()
                pdb = open(self.dirname + file, "w")
                pdb.write(read)
                pdb.close()
                input = str(self.dirname + file)
                PDBParser().get_structure(id, input)
                try:
                    test_structure = PDBParser().get_structure(id, input)
                    test_model = test_structure[0]
                    test_model[chain]
                except: 
                    raise StandardError, "%s - Invalid chain" % (chain)
            except:                             
                raise StandardError, "%s - Invalid PDB ID or PDB file" % (id)                            
        
    def sequencePDB(self, file, id, chain):
        "Extracts a sequence from the ATOM lines of a PDB file"
        
        # sequence from atom lines
        input = str(self.dirname + file)
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        string = ""
        for line in structure:
            if line[0:4] == "ATOM":
                if line[21] == str(chain):
                    CA = line[13:16]
                    res = line[17:20]
                    if CA == "CA ":
                        if res in aa_list:
                            string += aa_symbols[res]
                        else: pass
        sequence = string
                        
        output = str(self.dirname + id + ".fasta")
        out = open(output, "w")
        print >> out, ">Query_id" + "\n" + sequence + "\n"
        out.close()
        
        # full sequence from acession number on DBREF lines
        input = str(self.dirname + file)
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        
        for line in structure:
            if line[0:5] == "DBREF":
                if line[21] == str(chain):
                    data = line.split()
                    ch = data[2]
                    if ch == chain:
                        ac_number = data[6]
        try:
            fetch = Entrez.efetch(db="protein", id=ac_number, rettype="fasta")
                          
            output = str(self.dirname + id + ".fa")
            out = open(output, "w")
            out.write(fetch.read())
            out.close()
            read = SeqIO.parse(output)
            for record in read:
                sequence = str(record.seq)
            out = open(output, "w")
            print >> out, ">Query_id" + "\n" + sequence + "\n"
            out.close()
        except:
            copyfile(self.dirname + id + ".fasta", self.dirname + id + ".fa")

    def surfacePDB(self, file, id, chain):
        """"
        Points out surface residues in a PDB file (ASA > 7% (A^2))*
        *De et al.,2005. http://www.biomedcentral.com/1472-6807/5/15
        """
        
        input = str(self.dirname + file)
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        
        input_final = str("./" + file)
        out = open(input_final, "w")
        
        for line in structure:
            if line[0:4] == "ATOM":
                if line[21] == str(chain):
                    res = line[17:20]
                    res = res.rstrip()
                    res = res.lstrip()
                    if str(res) in aa_list:
                        print >> out, line.rstrip("\n")
        out.close()
        
        output = str("./" + file + ".txt")
        sasa.SASA(input_final, output)
        
        list = []
        input = str("./" + file + ".txt")
        op = open(input)
        read = op.readlines()
        for line in read:
            line = line.rstrip()
            line = line.split()
            if line[0] == chain:
                amino = str(line[1])
                res = int(line[2])
                area = float(line[3])
                if amino in aa_list:
                    info = [amino, res, area]
                    list.append(info)
        
        
        threshold = LP(self.parameterfile, "surface_threshold")
        surface = []
        asa_list = []
        total = 0
        for i in range(0, (len(list) - 1), 1):
            if list[i][0] == list[i + 1][0]:
                total += list[i][2]
            else:
                amino = str(list[i][0])
                res = list[i][1]
                area = total + list[i][2]
                value = [amino, res, area]
                asa_list.append(area)
                surface.append(value)
                total = 0
                pass
        
        output = str(self.dirname + id + ".surface")
        out = open(output, "w")
        
        asa_max = int(round(float(max(asa_list))))
        thrd = threshold * asa_max * 1.0 / 100
        for i in range(len(surface)):
            amino = str(surface[i][0])
            res = str(surface[i][1])
            area = float(surface[i][2])
            if area > thrd:
                print >> out, amino, res + "\t" + str(area)
        out.close()
        
        time.sleep(2)        
        try:
            remove("./" + file)
        except: pass
        try:
            remove("./" + file + ".txt")
        except: pass
        
    def parseSurfacePDB(self, id):
        "Parses residues at the surface level" 
        
        input = str(self.dirname + id + ".surface")
        input_surface = open(input, "r")
        surface = input_surface.readlines()
        input_surface.close()
        
        surface_points = []
        for line in surface:
            l = line.split()
            res_nb = int(l[1])
            surface_points.append(res_nb)
            
        return surface_points
        
    
    def matchResiduePosition(self, id, chain):
        "Gets residue positions for use in coevolution analysis"
        
        input = str(self.dirname + id + ".pdb")
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
                
        protein = []
        for line in structure:
            if line[0:4] == "ATOM":
                if line[21] == str(chain):
                    CA = line[13:16]
                    res_nb = line[22:26]
                    if CA == "CA ":
                        res_nb = line[22:26]
                        res = line[17:20]
                        res = res.rstrip()
                        res = res.lstrip()
                        if str(res) in aa_list:
                            protein.append(int(res_nb))
        return protein
    
    def copySequence(self, id):
        "Doubles the sequence files"
        
        copyfile(self.dirname + id + ".fa", self.dirname + id + "_1.fa")
        copyfile(self.dirname + id + ".fa", self.dirname + id + "_2.fa")
        remove(self.dirname + id + ".fa")
        return

