###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

import os
from Parameters import LoadParameters as LP
from src.UTILS import charge, charge_his, polarity, hydropathy
from os import remove, system
from shutil import copyfile
from itertools import combinations
from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline

class alignment:
    """
    Main code for multiple sequence alignment and scoring.
    
    Methods for computing MSAs:
    Clustalw - Chenna et al, 2003
    Muscle - Edgar, 2004
    Mafft - Katoh et al, 2002
    
    Methods for scoring MSAs:
    Sum-of-Pairs - Murata et al, 1985
    TODO: Circular Sum - Gonnet et al, 2000
    """
    def __init__(self, id1, id2, alignment, parameterfile, dirname):
        self.id1 = id1
        self.id2 = id2
        self.alignment = alignment
        self.parameterfile = parameterfile
        self.dirname = dirname
        
    def __call__(self, id1, id2, alignment, parameterfile, dirname):
        self.id1 = id1
        self.id2 = id2
        self.alignment = alignment
        self.parameterfile = parameterfile
        self.dirname = dirname
        

    def computeAlignment(self, id, alignment):
        "Computes multiple sequence alignment with inputed method"
        
        if alignment == "clustalw":
            gop = LP(self.parameterfile, "clustalw_gap_opening")
            gep = LP(self.parameterfile, "clustalw_gap_extension")
            d_matrix = LP(self.parameterfile, "clustalw_distance_matrix")
            
            input_sequences = self.dirname + id + ".fasta"
            output_align = self.dirname + id + ".aln"
            output_fasta = self.dirname + id + "_clustalw.fasta"
            output_tree = self.dirname + id + ".dnd"
            try:
                cmd = str(os.getcwd() + "/src/clustalw/clustalw.exe")
                clustalw = ClustalwCommandline(cmd, infile=input_sequences,
                                               outfile=output_align,
                                               newtree=output_tree,
                                               align="input",
                                               seqnos="ON",
                                               outorder="input",
                                               type="PROTEIN",
                                               pwmatrix=d_matrix,
                                               gapopen=gop,
                                               gapext=gep) 
                clustalw()
            except:
                cmd = str(os.getcwd() + "/src/clustalw/clustalw")
                clustalw = ClustalwCommandline(cmd, infile=input_sequences,
                                               outfile=output_align,
                                               newtree=output_tree,
                                               align="input",
                                               seqnos="ON",
                                               outorder="input",
                                               type="PROTEIN",
                                               pwmatrix=d_matrix,
                                               gapopen=gop,
                                               gapext=gep) 
                clustalw()
            AlignIO.convert(output_align, "clustal", output_fasta, "fasta")
            try:
                remove(output_align)
                remove(output_tree)
            except:
                pass
            
        elif alignment == "muscle":
            iteration = LP(self.parameterfile, "muscle_max_iteration")
            
            input_sequences = self.dirname + id + ".fasta"
            output_align = self.dirname + id + "_muscle.aln"
            output_fasta = self.dirname + id + "_muscle.fasta"
            
            muscle = MuscleCommandline(input=input_sequences,
                                       out=output_align,
                                       clwstrict=True,
                                       maxiters=iteration)
            muscle()
            AlignIO.convert(output_align, "clustal", output_fasta, "fasta")
            try:
                remove(output_align)
            except:
                pass
            
            organism_order = []
            input_sequences = self.dirname + id + ".fasta"
            align = SeqIO.parse(input_sequences, "fasta", IUPAC.protein)
            for record in align:
                org = record.description
                organism_order.append(org)
                
            rec = dict()
            output_fasta = self.dirname + id + "_muscle.fasta"
            align = SeqIO.parse(output_fasta, "fasta", IUPAC.protein)
            for record in align:
                org = str(record.description)
                seq = str(record.seq)
                rec[org] = seq
            
            fasta = open(output_fasta, "w")
            fasta.close()
            fasta = open(output_fasta, "a")
            for org in (organism_order):
                seq = rec[org]
                fasta.write(">" + org + "\n" + seq + "\n")
            fasta.close()
            
        else:
            configuration = LP(self.parameterfile, "mafft_configuration")
            threads = LP(self.parameterfile, "mafft_threading")
            input_sequences = self.dirname + id + ".fasta"
            output_fasta = self.dirname + id + "_mafft.fasta"
            
            if configuration == "fftnsi":
                if threads == False:
                    fftnsi = "mafft --retree 2 --maxiterate 1000 --inputorder "
                    mafft = system(fftnsi + input_sequences + ">" + output_fasta)
                    mafft
                else:
                    try:
                        threads = int(threads)
                        fftnsi = "mafft --retree 2 --maxiterate 1000\
                         --inputorder --threads %i " % (threads)
                        mafft = system(fftnsi + input_sequences + ">" + output_fasta)
                        mafft
                    except:
                        fftnsi = "mafft --retree 2 --maxiterate 1000 --inputorder "
                        mafft = system(fftnsi + input_sequences + ">" + output_fasta)
                        mafft
            else:
                if threads == False:
                    linsi = "mafft --localpair --maxiterate 1000 --inputorder "
                    mafft = system(linsi + input_sequences + ">" + output_fasta)
                    mafft
                else:
                    try:
                        threads = int(threads)
                        linsi = "mafft --localpair --maxiterate 1000\
                         --inputorder --threads %i " % (threads)
                        mafft = system(linsi + input_sequences + ">" + output_fasta)
                        mafft
                    except:
                        linsi = "mafft --localpair --maxiterate 1000 --inputorder "
                        mafft = system(linsi + input_sequences + ">" + output_fasta)
                        mafft    
        
    def cutAlignment(self, file, id, alignment):
        "Selects MSA columns of interest (Query_id != '-')"
               
        description = []
        align = []
        columns = []
        positions = []
        blocks = []
        new_align = []
        new_align_ord = []
        new_align_concate = []
        self.cut_alignment = []
        aa_red = LP(self.parameterfile, "alphabet_reduction")
        
        if alignment != "custom":
            input = self.dirname + id + "_" + alignment + ".fasta"
            alignment = AlignIO.read(input, "fasta")
            for record in alignment:
                key = record.id
                description.append(key)
            
            k = int(-1)
            for s in description:
                k += 1
                key = s.find("Query_id")
                if key != -1:
                    break
            
            align_length = alignment.get_alignment_length()
            for position in range(0, align_length):
                column = alignment[:, position]
                align.append(column)
                if column[k] != "-":
                    columns.append(column)
                    positions.append(position)
                    
            for i in range(0, len(positions), 1):
                beg = int(positions[i])
                end = int(positions[i] + 1)
                block = alignment[:, beg:end]         
                blocks.append(block)
            
            for block in blocks:
                for record in block:
                    seq = str(record.seq)
                    new_align.append(seq)
            
            numb_blocks = len(new_align) / len(columns[0])
            for i in range(0, len(columns[0])):
                for j in range(0, len(new_align), len(columns[0])):
                    new_align_ord.append(new_align[i + j])
        
            for i in range(0, len(new_align_ord), numb_blocks):
                pseudolist = new_align_ord[i:i + numb_blocks]
                list = ""
                for j in pseudolist:
                    list += j
                new_align_concate.append(list)
    
            for seq in new_align_concate:
                if aa_red != False:
                    red = [AR(e, aa_red) for e in seq]
                    new_seq = ""
                    for i in red:
                        new_seq += str(i)
                    self.cut_alignment.append(new_seq)
                else:
                    self.cut_alignment.append(seq) 
                
            return self.cut_alignment
        
        else:
            output = self.dirname + id + "_" + alignment + ".fasta"
            copyfile(self.dirname + file, output)
            input = self.dirname + id + "_" + alignment + ".fasta"
           
            alignment = AlignIO.read(input, "fasta")
            for record in alignment:
                key = record.id
                description.append(key)
            
            align_length = alignment.get_alignment_length()
            for position in range(0, align_length):
                column = alignment[:, position]
                align.append(column)
                if column[0] != "-":
                    columns.append(column)
                    positions.append(position)
                    
            for i in range(0, len(positions), 1):
                beg = int(positions[i])
                end = int(positions[i] + 1)
                block = alignment[:, beg:end]         
                blocks.append(block)
            
            for block in blocks:
                for record in block:
                    seq = str(record.seq)
                    new_align.append(seq)
            
            numb_blocks = len(new_align) / len(columns[0])
            for i in range(0, len(columns[0])):
                for j in range(0, len(new_align), len(columns[0])):
                    new_align_ord.append(new_align[i + j])
        
            for i in range(0, len(new_align_ord), numb_blocks):
                pseudolist = new_align_ord[i:i + numb_blocks]
                list = ""
                for j in pseudolist:
                    list += j
                new_align_concate.append(list)
    
            for seq in new_align_concate:
                if aa_red != False:
                    red = [AR(e, aa_red) for e in seq]
                    new_seq = ""
                    for i in red:
                        new_seq += str(i)
                    self.cut_alignment.append(new_seq)
                else:
                    self.cut_alignment.append(seq)      
               
            return self.cut_alignment
            
     
    def alignScore(self, id, alignment):
        """
        Computes a score for the MSA inputed.
        
        Methods implemented:
        Sum-of-pairs (SP) score -  Murata et al, 1985
        as explained in Gonnet et al, 2000. 
        SP is the sum of all possible combinations of
        pairwise scores. 
        
        !!Disclaimer: alignmentScore is terribly slow!!
        
        (To Do - Circular Sum by Gonnet et al, 2000)
        """
        score = LP(self.parameterfile, "alignment_score")
        
        if score == "sumofpairs":
            input = self.dirname + id + "_" + alignment + ".fasta"
            sequences = []
            input_sequences = SeqIO.parse(input, "fasta", IUPAC.protein)
            for record in input_sequences:
                seq = str(record.seq)
                sequences.append(seq) 
            
            SumOfPairs = 0
            for pair in combinations(sequences, 2): 
                SumOfPairs += pairwiseScore(pair[0], pair[1])
            
            print SumOfPairs
        else: pass

def pairwiseScore(seq1, seq2):
    """
    s(x,y) = { matchScore(x,y) if x!='-' and y!='-';
               0 if x=='-' and y=='-', because the delection as caused earlier
               gap penalty, depending on the gap length  gap + length * increment}
               
    gap - depends on the scoring matrix (PAM, BLOSUM, etc)
    length - length of the gap
    increment - incremental penalty that depends on the scoring matrix
    
    BLOSUM62, gap = -4, increment = 1 -> increment = length 
    """
    
    gap = -4.0
    incr_top = 0
    incr_bottom = 0
    pairwise_score = 0
    for i, j in zip(range(len(seq1)), range(len(seq2))):
        aa1 = seq1[i]
        aa2 = seq2[j] 
        if aa1 == "-" and aa2 == "-" :
            pairwise_score += 0
        elif aa1 != "-" and aa2 != "-":
            pairwise_score += float(matchScore(aa1, aa2, "BLOSUM62"))
        elif aa1 == "-" and aa2 != "-":
            try:
                aa11 = seq1[i + 1]
                aa22 = seq2[j + 1]
                if aa11 == "-" and aa22 != "-":
                    incr_top += 1
                else: 
                    pairwise_score += gap + incr_top * incr_top
                    incr_top = 0
            except: 
                pairwise_score += gap
                pass
        elif aa1 != "-" and aa2 == "-":
            try:
                aa11 = seq1[i + 1]
                aa22 = seq2[j + 1]
                if aa11 != "-" and aa22 == "-":
                    incr_bottom += 1
                else: 
                    pairwise_score += gap + incr_bottom * incr_bottom
                    incr_bottom = 0
            except: 
                pairwise_score += gap
                pass
        else: pass
        
    return pairwise_score
         
def matchScore(alpha, beta, matrix):
    "Matches scores from a matrix"
    
    alphabet = {}    
    alphabet["A"] = 0
    alphabet["R"] = 1
    alphabet["N"] = 2
    alphabet["D"] = 3
    alphabet["C"] = 4
    alphabet["Q"] = 5
    alphabet["E"] = 6
    alphabet["G"] = 7
    alphabet["H"] = 8
    alphabet["I"] = 9
    alphabet["L"] = 10
    alphabet["K"] = 11
    alphabet["M"] = 12
    alphabet["F"] = 13
    alphabet["P"] = 14
    alphabet["S"] = 15
    alphabet["T"] = 16
    alphabet["W"] = 17
    alphabet["Y"] = 18
    alphabet["V"] = 19
    alphabet["B"] = 20
    alphabet["Z"] = 21
    alphabet["X"] = 22
    alphabet["-"] = 22
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return mapMatrix(matrix)[lut_x][lut_y]
    
def mapMatrix(matrix):
    "Maps a matrix of floats"
    matrix = matrix.upper()
    
    score_matrix = []
    input = './Matrix/' + matrix
    input_matrix = open(input, 'r')
    for line in input_matrix.readlines():
        score_matrix.append(map(float, line.split()))
    input_matrix.close()
    
    return score_matrix
    
def AR(aminoacid, method):
    """Performs alphabet reduction.
    Alphabets: charge, charge_his, polarity, hydropathy
    """

    if method == "charge":
        return charge[aminoacid]
    elif method == "charge_his":
        return charge_his[aminoacid]
    elif method == "polarity":
        return polarity[aminoacid]
    elif method == "hydropathy":
        return hydropathy[aminoacid]
    else:
        return aminoacid
    
        
