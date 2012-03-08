###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from src.SEQ import sequence
from src.BLAST import psiblast
from src.ORGANISM import organism
from src.ALIGN import alignment
from src.COEVOL import coevolution
from src.INFO import information
from Parameters import results_sifts

class main:
    """
    Main script caller
    """
    def __init__(self, file1, file2, id1, id2, chain1, chain2, 
                 psiblast, alignment, coevolution):
        self.file1 = file1
        self.file2 = file2
        self.id1 = id1
        self.id2 = id2
        self.chain1 = chain1
        self.chain2 = chain2
        self.psiblast = psiblast
        self.alignment = alignment
        self.coevolution = coevolution
    def __call__(self, file1, file2, id1, id2, chain1, chain2, 
                 psiblast, alignment, coevolution):
        self.file1 = file1
        self.file2 = file2
        self.id1 = id1
        self.id2 = id2
        self.chain1 = chain1
        self.chain2 = chain2
        self.psiblast = psiblast
        self.alignment = alignment
        self.coevolution = coevolution
    
    
    def sequenceSripts(self):
        seq = sequence(self.file1, self.file2, self.id1, self.id2, 
                       self.chain1, self.chain2)
        if self.id1 != self.id2:
            if self.chain1 == "" and self.chain2 == "":
                seq.validFASTA(self.file1, self.id1)
                seq.queryFASTA(self.file1, self.id1)
                seq.validFASTA(self.file2, self.id2)
                seq.queryFASTA(self.file2, self.id2)
            else:
                seq.validPDB(self.file1, self.id1, self.chain1)
                seq.sequencePDB(self.file1, self.id1, self.chain1)
                seq.surfacePDB(self.file1, self.id1, self.chain1)
                seq.validPDB(self.file2, self.id2, self.chain2)
                seq.sequencePDB(self.file2, self.id2, self.chain2)
                seq.surfacePDB(self.file2, self.id2, self.chain2)
        else:
            if self.chain1 == "" and self.chain2 == "":
                seq.validFASTA(self.file1, self.id1)
                seq.queryFASTA(self.file1, self.id1)
            else:
                if self.chain1 != self.chain2:
                    seq.validPDB(self.file1, self.id1, self.chain1)
                    seq.sequencePDB(self.file1, self.id1 + "_1", self.chain1)
                    seq.surfacePDB(self.file1, self.id1 + "_1", self.chain1)
                    seq.validPDB(self.file1, self.id1, self.chain2)
                    seq.sequencePDB(self.file1, self.id1 + "_2", self.chain2)
                    seq.surfacePDB(self.file1, self.id1 + "_2", self.chain2)
                else:
                    seq.validPDB(self.file1, self.id1, self.chain1)
                    seq.sequencePDB(self.file1, self.id1, self.chain1)
                    seq.surfacePDB(self.file1, self.id1, self.chain1)
        return
    
    def psiblastSripts(self):
        seq = sequence(self.file1, self.file2, self.id1, self.id2, 
                       self.chain1, self.chain2)
        blast = psiblast(self.id1, self.id2, self.psiblast)
        if self.id1 != self.id2:
            blast.searchPSIBLAST(self.id1,self.psiblast)
            blast.searchPSIBLAST(self.id2,self.psiblast)
            blast.validXML(self.id1)
            blast.validXML(self.id2)
            blast.sequencesXML(self.id1,self.psiblast)
            blast.sequencesXML(self.id2,self.psiblast)
        else:
            if self.chain1 == "" and self.chain2 == "":
                seq.copySequence(self.id1)
                blast.searchPSIBLAST(self.id1 + "_1",self.psiblast)
                blast.searchPSIBLAST(self.id1 + "_2",self.psiblast)
                blast.validXML(self.id1 + "_1")
                blast.validXML(self.id1 + "_2")
                blast.sequencesXML(self.id1 + "_1",self.psiblast)
                blast.sequencesXML(self.id1 + "_2",self.psiblast)
            else:
                if self.chain1 != self.chain2:
                    blast.searchPSIBLAST(self.id1 + "_1",self.psiblast)
                    blast.searchPSIBLAST(self.id1 + "_2",self.psiblast)
                    blast.validXML(self.id1 + "_1")
                    blast.validXML(self.id1 + "_2")
                    blast.sequencesXML(self.id1 + "_1",self.psiblast)
                    blast.sequencesXML(self.id1 + "_2",self.psiblast)
                else:
                    seq.copySequence(self.id1)
                    blast.searchPSIBLAST(self.id1 + "_1",self.psiblast)
                    blast.searchPSIBLAST(self.id1 + "_2",self.psiblast)
                    blast.validXML(self.id1 + "_1")
                    blast.validXML(self.id1 + "_2")
                    blast.sequencesXML(self.id1 + "_1",self.psiblast)
                    blast.sequencesXML(self.id1 + "_2",self.psiblast)
        return
    
    def organismSripts(self):
        org = organism(self.id1, self.id2, self.psiblast)
        if self.id1 != self.id2:
            org.uniqueOrganism(self.id1, self.id2)
            org.pairwiseDistance(self.id1, self.id2)
            org.getsCorrelation()
            org.removeSequences(self.id1, self.id2)
        else:
            org.uniqueOrganism(self.id1 + "_1", self.id1 + "_2")
            org.pairwiseDistance(self.id1 + "_1", self.id1 + "_2")
            org.getsCorrelation()
            org.removeSequences(self.id1 + "_1", self.id1 + "_2")                               
        return
    
    def alignmentSripts(self): 
        aln = alignment(self.id1, self.id2, self.alignment)
        if self.id1 != self.id2:
            aln.computeAlignment(self.id1, self.alignment)
            aln.computeAlignment(self.id2, self.alignment)
            #aln.alignScore(self.id1, self.alignment)
            #aln.alignScore(self.id2, self.alignment)
        else:
            aln.computeAlignment(self.id1 + "_1", self.alignment)
            aln.computeAlignment(self.id1 + "_2", self.alignment)
            #aln.alignScore(self.id1 + "_1", self.alignment)
            #aln.alignScore(self.id1 + "_2", self.alignment)                              
        return
    
    def coevolutionSripts(self):
        coevol = coevolution(self.file1, self.file2, self.id1, self.id2, 
                             self.chain1, self.chain2, 
                             self.alignment, self.coevolution)
        if self.id1 != self.id2:
            coevol.coevolAnalysis(self.file1, self.file2,
                                  self.id1, self.id2, 
                                  self.chain1, self.chain2, 
                                  self.alignment, self.coevolution)
            coevol.bestInfo(self.id1, self.id2,  
                                  self.alignment, self.coevolution)
            if self.chain1 == "" and self.chain2 == "":
                pass
            else:
                coevol.structureSingle(self.id1, self.id2, 
                                            self.chain1, self.chain2, 
                                            self.alignment, self.coevolution)
                
        else:
            coevol.coevolAnalysis(self.file1, self.file1,
                                  self.id1 + "_1", self.id1 + "_2", 
                                  self.chain1, self.chain2, 
                                  self.alignment, self.coevolution)
            coevol.bestInfo(self.id1 + "_1", self.id1 + "_2",  
                                  self.alignment, self.coevolution)
            if self.chain1 == "" and self.chain2 == "":
                pass
            else:
                if self.chain1 != self.chain2:
                    coevol.structurePair(self.id1, self.id1, 
                                            self.chain1, self.chain2, 
                                            self.alignment, self.coevolution)
        return
    
    def infoScripts(self):
        info = information(self.id1, self.id2,self.chain1, self.chain2)
        if self.id1 != self.id2:
            if self.chain1 == "" and self.chain2 == "":
                info.getInfo(self.id1)
                info.getInfo(self.id2)
            else:
                info.getInfo(self.id1)
                info.getInfo(self.id2)
                if results_sifts == True:
                    info.getSIFTS(self.id1, self.chain1)
                    info.getSIFTS(self.id2, self.chain2)
                else: pass
        else:
            if self.chain1 == "" and self.chain2 == "":
                info.getInfo(self.id1 + "_1")
            else:
                if self.chain1 != self.chain2:
                    info.getInfo(self.id1 + "_1")
                    info.getInfo(self.id1 + "_2")
                    if results_sifts == True:
                        info.getSIFTS(self.id1 + "_1", self.chain1)
                        info.getSIFTS(self.id1 + "_2", self.chain2)
                    else: pass
                else:
                    info.getInfo(self.id1 + "_1")
                    if results_sifts == True:
                        info.getSIFTS(self.id1 + "_1", self.chain1)
                    else: pass
        return

    
    
