###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from Bio import SeqIO

class information:
    """
    Main code for generating extended results.
    """
    def __init__(self, id1, id2, chain1, chain2, dirname):
        self.id1 = id1
        self.id2 = id2
        self.chain1 = chain1
        self.chain2 = chain2
        self.dirname = dirname
                
    def __call__(self, id1, id2, chain1, chain2, dirname):
        self.id1 = id1
        self.id2 = id2
        self.chain1 = chain1
        self.chain2 = chain2
        self.dirname = dirname
    
    def getInfo(self, id):
        "Creates info about the sequences, psiblast, organisms, etc"
    
        input = self.dirname + id + ".fasta"
        sequence = SeqIO.parse(input, "fasta")
        for seq_record in sequence:
            seq = seq_record.seq
            length = len(seq)
            break
    
        input = self.dirname + id + ".blast"
        hit = 0
        sequences = SeqIO.parse(input, "fasta")
        for record in sequences:
            sequence = record.seq              
            hit += 1      
        
        input = self.dirname + id + ".fasta"
        sequences = SeqIO.parse(input, "fasta")
        organisms = 0
        for record in sequences:
            sequence = record.seq
            organisms += 1

        
        output = self.dirname + "results.txt"
        out = open(output, "a")
        print >> out, "ID" + "\t" + "LengSeq" + "\t" + "NHits" + "\t" + \
        "NOrganisms"
        print >> out, str(id) + "\t" + str(length) + "\t" + \
        str(hit) + "\t" + str(organisms) + "\n"
        out.close()
        
    def getSIFTS(self, id, chain):
        """
        Web_Services based on SIFTS @ 
        http://www.ebi.ac.uk/pdbe/docs/sifts/
        """
        
        id = id.lower()
        try:
            id = id.rstrip("_1")
        except:
            pass
        try:
            id = id.rstrip("_2")
        except:
            pass
        
        # Uniprot ID and SCOP
        input = "./SIFTS/pdb_chain_scop_uniprot.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()
        
        unip = "Not_found"
        scop = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    unip = str(l[2])
                    scop = str(l[5])
        
        # CATH
        input = "./SIFTS/pdb_chain_cath_uniprot.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()

        cath = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    cath = str(l[4])

        
        # EC (enzyme)
        input = "./SIFTS/pdb_chain_enzyme.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()
        
        enz = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    enz = str(l[4])
        
        # Interpro
        input = "./SIFTS/pdb_chain_interpro.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()
        
        inter = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    inter = str(l[2])
        
        # Pfam
        input = "./SIFTS/pdb_chain_pfam.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()
        
        pfam = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    pfam = str(l[4])
 
        # Taxonomy
        input = "./SIFTS/pdb_chain_taxonomy.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()

        taxid = "Not_found"
        taxnm = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                if l[1] == str(chain):
                    taxid = str(l[2])
                    taxnm = str(l[7])

        
        # Pubmed
        input = "./SIFTS/pdb_pubmed.lst"
        sifts = open(input, "r")
        read = sifts.readlines()
        sifts.close()

        pubm = "Not_found"
        for line in read:
            if line[0:4] == str(id):
                l = line.rstrip("\n")
                l = l.split("\t")
                pubm = str(l[2])
        
        output = self.dirname + "bioresults.txt"
        out = open(output, "a")
        print >> out, "Protein_ID" + "\t" + "Uniprot" + "\t" + "SCOP" + "\t" + \
        "CATH" + "\t" + "Enzyme_EC" + "\t" + "Interpro" + "\t" + "Pfam" + "\t" + \
        "Taxonomy_id" + "\t" + "Taxonomy_name" + "\t" + "Pubmed" 
        
        print >> out, str(id) + "\t" + str(unip) + "\t" + \
        str(scop) + "\t" + str(cath) + "\t" + str(enz) + "\t" + \
        str(inter) + "\t" + str(pfam) + "\t" + str(taxid) + "\t" + \
        str(taxnm) + "\t" + str(pubm) + "\n"
        out.close()
        
        
