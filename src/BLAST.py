###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from Parameters import psiblast_evalue, psiblast_identity, psiblast_coverage
from Bio import SeqIO, Entrez
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Blast.Applications import NcbipsiblastCommandline
Entrez.email = "entrez@mail.com"


class psiblast:
    """
    Main code for psiblast search over internet or at local database.
    
    Method for searching homologous sequences:
    PSI-Blast - Altschul et al, 1997
    """
    def __init__(self, id1, id2, psiblast):
        self.id1 = str(id1)
        self.id2 = str(id2)
        self.psiblast = psiblast
    def __call__(self, id1, id2, psiblast):
        self.id1 = str(id1)
        self.id2 = str(id2)
        self.psiblast = psiblast

    def searchPSIBLAST(self, id, psiblast):
        "Psi-Blast over a local database or over the internet"
        
        if psiblast == "local":
            # edit psiblast_evalue at Parameters.py
            evalue = psiblast_evalue
            reference_protein = "refseq_protein"
        
            sequence = "./Data/" + id + ".fasta"
            output = "./Data/"+ id + ".xml" 
            psiblast = NcbipsiblastCommandline(query=sequence, 
										 db=reference_protein, 
										 outfmt=5,  
										 threshold=evalue, 
										 out=output) 
            psiblast()
        else:
            # edit psiblast_evalue at Parameters.py
            evalue = psiblast_evalue
            reference_protein = "refseq_protein"
        
            for seq_record in SeqIO.parse("./Data/" + id + ".fasta", 
                                          "fasta",IUPAC.protein):
                sequence = seq_record.seq
        
                psiblast = NCBIWWW.qblast("blastp", 
								    reference_protein, 
								    sequence, 
								    service="psi", 
								    expect=evalue,
								    hitlist_size=500)
                psiblast

            output = "./Data/"+ id + ".xml"
            saveblast = open(output, "w")
            saveblast.write(psiblast.read())
            saveblast.close()
            psiblast.close()

    def validXML(self, id):
        "Checks if the input file is a valid XML"
        
        try:
            input ="./Data/" + id + ".xml"
            input_xml = open(input, "r")
            xml = input_xml.readline()
            input_xml.close()
            if xml[0:5] == "<?xml":
                pass
            else:                             
                raise StandardError, "%s - Invalid xml" %(input)
        except:
            raise StandardError, "%s - Invalid xml or not found" %(input)

    def sequencesXML(self, id):
        "Extracts records from xml and writes FASTA (full-length) sequences"
        
        # edit psiblast_identity and psiblast_coverage at Parameters.py
        thresh_identity = psiblast_identity
        thresh_coverage = psiblast_coverage
        
        input ="./Data/" + id + ".xml"
        input_xml = open(input, "r")
        
        hits = []
        for record in NCBIXML.parse(input_xml):            
            for align in record.alignments:
                hit_id = align.hit_id
                for hsp in align.hsps:    
                    positives = int(hsp.positives)
                    identities = int(hsp.identities)
                    q_start = int(hsp.query_start)
                    q_end = int(hsp.query_end)
                    query = (q_end - q_start) * 1.0
                    sbjct1 = positives * 1.0
                    coverage = sbjct1 / query * 100
                    sbjct2 = identities * 1.0
                    identity =sbjct2 / query * 100
                    if coverage > thresh_coverage and identity > thresh_identity:
                        hits.append(hit_id)
        input_xml.close()
        
        if hits == []:
            raise StandardError, "%s - No Hits found in PSI-BLAST search" %(input) 
                 
        for hit_id in hits:
            gi = hit_id[hit_id.find("id|") + 4:hit_id.find("|ref")]
            efetch = Entrez.efetch(db="protein", id=gi, rettype="fasta")
            for values in efetch:
                description = values
                break
            sequence = ""
            for values in efetch:
                sequence += values.rstrip("\n")
            try: 
                organism = description[description.find("[") + 1:description.find("]")]
                organism = organism.split()
                if len(organism) != 1:
                    species = str(organism[0] + "_" + organism[1])
                else:
                    species = str(organism[0] + "_" + "sp.")
                output = "./Data/" + id + ".blast"
                blast = open(output, "a")
                blast.write("\n" + ">" + species + "\n" + sequence + "\n")
                blast.close()
            except: 
                raise StandardError, "%s - No Hits found in PSI-BLAST search" %(input)    
            
 
        
