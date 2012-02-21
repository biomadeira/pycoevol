###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from src.UTILS import aa_list, aa_symbols
from Parameters import surface_threshold, interface_threshold
from os import system, remove
from shutil import copyfile
from Bio import SeqIO, Entrez
from Bio.Alphabet import IUPAC
from Bio.PDB.PDBParser import PDBParser
Entrez.email = "entrez@mail.com"


class sequence:
    """
    Main code for handling sequences and structures.
    """
    def __init__(self, file1, file2, id1, id2, chain1, chain2):
        self.file1 = str(file1)
        self.file2 = str(file2)
        self.chain1 = str(chain1)
        self.chain2 = str(chain2)
        self.id1 = str(id1)
        self.id2 = str(id2)
    def __call__(self, file1, file2, id1, id2, chain1, chain2):
        self.file1 = str(file1)
        self.file2 = str(file2)
        self.chain1 = str(chain1)
        self.chain2 = str(chain2)
        self.id1 = str(id1)
        self.id2 = str(id2)
       
    def validFASTA(self, file):
        "Checks if the input file is a valid FASTA file"
        
        try:
            input = str("./Data/" + file)
            SeqIO.read(input, "fasta", IUPAC.protein)
        except:                           
            raise StandardError, "%s - Invalid sequence" %(file)
    
    def queryFASTA(self, file):
        "Changes FASTA original header to 'Query_id'"
        
        input = str("./Data/" + file)
        input_sequence = SeqIO.parse(input, "fasta", IUPAC.protein)
        for record in input_sequence:
            sequence = str(record.seq)
            break
            
        output = str("./Data/" + id + ".fa")
        out = open(output, "w")
        print >> out, ">Query_id" + "\n" + sequence + "\n"
        out.close()
        
        remove(input)
        
    def validPDB(self, file, id, chain):
        "Checks if the input file is a valid PDB file"
        
        try:
            input = str("./Data/" + file)
            PDBParser().get_structure(id, input)
            try:
                test_structure = PDBParser().get_structure(id, input)
                test_model = test_structure[0]
                test_model[chain]
            except: 
                raise StandardError, "%s - Invalid chain" %(chain)
        except:                             
            raise StandardError, "%s - Invalid structure" %(file)
        
    def sequencePDB(self, file, id, chain):
        "Extracts a sequence from the ATOM lines of a PDB file"
        
        # sequence from atom lines
        input = str("./Data/" + file)
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        string=""
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
                        
        output = str("./Data/" + id + ".fasta")
        out = open(output, "w")
        print >> out, ">Query_id" + "\n" + sequence + "\n"
        out.close()
        
        # full sequence from acession number on DBREF lines
        input = str("./Data/" + file)
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
                          
            output = str("./Data/" + id + ".fa")
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
            copyfile("./Data/" + id + ".fasta", "./Data/" + id + ".fa")

    def surfacePDB(self, file, id, chain):
        """"
        Points out surface residues in a PDB file (ASA > 7% (A^2))*
        *De et al.,2005. http://www.biomedcentral.com/1472-6807/5/15
        """
        
        input = str("./Data/" + file)
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        
        output = str("./" + file)
        out = open(output, "w")
        
        for line in structure:
            if line[0:4] ==  "ATOM":
                if line[21] == str(chain):
                    res = line[17:20]
                    res = res.rstrip()
                    res = res.lstrip()
                    if str(res) in aa_list:
                        print >> out, line.rstrip("\n")
        out.close()

        pdbsurface = system("pdbsurface")
        pdbsurface
        
        input = str("./" + file + ".txt")
        input_surface = open(input, "r")
        surface = input_surface.readlines()
        input_surface.close()
        
        output = str("./Data/" + id + ".surface")
        out = open(output, "w")
        
        # edit surface_threshold at Parameters.py
        threshold = surface_threshold
        asa_list = []
        for line in surface:
            l = line.split()
            asa = float(l[2])
            asa_list.append(asa)
        
        asa_max = max(asa_list)
        thrd = threshold * asa_max * 1.0 / 100
        for line in surface:
            l = line.split()
            asa = float(l[2])
            asa_list.append(asa)
            if asa  > thrd:
                print >> out, line.rstrip("\n")
        out.close()
                
        try:
            remove("./" + file)
            remove("./" + file + ".txt")
        except: pass
        

    def interfacePDB(self, file, id, chain1, chain2):
        """
        Points out interface contact points if file1==file2, chain1!=chain2, 
        and chains are in contact.
        Interface residues are defined as a point at the surface which lose 
        solvent-accessible surface area (by >1 A^2)* on oligomerisation.
        *De et al.,2005. http://www.biomedcentral.com/1472-6807/5/15
        """ 
            
        input = str("./Data/" + file)
        structure = PDBParser().get_structure(id, input)
        model = structure[0]
        chain_1 = model[chain1]
        chain_2 = model[chain2]
        
        input_structure = open(input, "r")
        structure = input_structure.readlines()
        input_structure.close()
        
        protein1 = []
        protein2 = []
        for line in structure:
            if line[0:4] ==  "ATOM":
                if line[21] == chain1:
                    if line[13:16] == "CA ":
                        res = line[17:20]
                        res = res.rstrip()
                        res = res.lstrip()
                        res_nb = line[22:26]
                        res_nb = res_nb.rstrip()
                        res_nb = res_nb.lstrip()
                        if str(res) in aa_list:
                            protein1.append(int(res_nb))
                if line[21] == chain2:
                    if line[13:16] == "CA ":
                        res = line[17:20]
                        res = res.rstrip()
                        res = res.lstrip()
                        res_nb = line[22:26]
                        res_nb = res_nb.rstrip()
                        res_nb = res_nb.lstrip()
                        if str(res) in aa_list:
                            protein2.append(int(res_nb))
         
        for res_one in protein1:
            for res_two in protein2:
                res1 = chain_1[res_one]
                res2 = chain_2[res_two]
                distance = residueDist(res1, res2)
                if distance < 9:
                    contacts = True
                    break
         
        if contacts == True:
            print "True"
            input = str("./Data/" + file)
            input_structure = open(input, "r")
            structure = input_structure.readlines()
            input_structure.close()
        
            output = str("./" + file)
            out = open(output, "w")
        
            for line in structure:
                if line[0:4] ==  "ATOM":
                    if line[21] == chain1 or line[21] == chain2:
                        res = line[17:20]
                        res = res.rstrip()
                        res = res.lstrip()
                        if str(res) in aa_list:
                            print >> out, line.rstrip("\n")
                elif line[0:6] == "HETATM":
                    pass
                else:
                    print >> out, line.rstrip("\n")
                           
            out.close()

            pdbsurface = system("pdbsurface -c")
            pdbsurface
        
            input = str("./" + id + ".contacts")
            input_interface = open(input, "r")
            interface = input_interface.readlines()
            input_interface.close()
        
            output = str("./Data/" + id + ".contacts")
            out = open(output, "w")
        
            # edit interface_threshold at Parameters.py
            threshold = interface_threshold
            for line in interface:
                l = line.split()
                asa = float(l[7])
                if asa  > threshold:
                    print >> out, line.rstrip("\n")
            out.close()
                
            try:
                remove("./" + file)
                remove("./" + id + ".contacts")
            except: pass
        else: pass
        
    def parseSurfacePDB(self, id):
        "Parses residues at the surface level" 
        
        input = str("./Data/" + id + ".surface")
        input_surface = open(input, "r")
        surface = input_surface.readlines()
        input_surface.close()
        
        surface_points = []
        for line in surface:
            l = line.split()
            res_nb = int(l[1])
            surface_points.append(res_nb)
            
        return surface_points
        
    def parseInterfacePDB(self, id):
        "Parses the interface contacts at the surface level" 
        
        input = str("./Data/" + id + ".contacts")
        input_interface = open(input, "r")
        interface = input_interface.readlines()
        input_interface.close()
        
        interface_points = []
        interface_1 = []
        interface_2 = []
        for line in interface:
            l = line.split()
            res_nb_1 = int(l[1])
            res_nb_2 = int(l[1])
            val = [res_nb_1, res_nb_2]
            interface_points.append(val)
            interface_1.append(res_nb_1)
            interface_2.append(res_nb_2)
            
        return interface_points, interface_1, interface_2
    
    def matchResiduePosition(self, id, chain):
        "Gets residue positions for use in coevolution analysis"
        
        input = str("./Data/" + id + ".pdb")
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
    
    def copySequence(self, id1):
        "Doubles the sequence files"
        
        copyfile("./Data/" + id + ".fasta", "./Data/" + id + "_1.fasta")
        copyfile("./Data/" + id + ".fasta", "./Data/" + id + "_2.fasta")
        remove("./Data/" + id + ".fasta")
        return
        
def residueDist(res_one, res_two): 
    "Returns the C-alpha distance between two residues"
    try: 
        resone = res_one["CA"]
        restwo = res_two["CA"]
        distance  = resone - restwo
        return distance
    except:  pass



    