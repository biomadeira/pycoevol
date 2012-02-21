###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

from src.SEQ import sequence as class_sequence
from src.ALIGN import alignment as class_alignment
from Parameters import results_histogram, results_heatmap, results_structure
from Parameters import best_results
from src.UTILS import aa
from math import log, e, factorial
from numpy import mean, std, zeros, sqrt
from matplotlib import pyplot


class coevolution:
    """
    Main code for coevolution analysis.
    
    Matrix-based Methods:
    * Residue Contact Preferences, Volume Normalized - Glaser et al, 2001.
    * Contact PDB-derived Likelihood Matrix  - Singer et al, 2002.
    * Residue-residue volume normalized  - based on Esque et al, 2010.
    
    Mutual Information based methods:
    * Mutual Information - Gloor el al, 2005.
    * MI by pair entropy - Martin el al, 2005.
    * Row and column weighed MI - Gouveia-Oliveira et al, 2007.
    * Contact preferences, volume normalized MIE - F. Madeira, 2012.
    (unpublished)
    
    Correlation-based methods:
    * OMES (Observed Minus Expected Squared) - Kass and Horovitz, 2002.
    * Pearson's correlation - Gobel et al, 1994. (slow)
    * Spearman's rank correlation - Pazos et al, 1997. (slow)
    * McBASC (McLachlan Based Substitution Correlation) - Fodor and 
    Aldrich, 2004. (slow)
    * Quartets - Galitsky, 2002. 
    
    Perturbation-based methods:
    * SCA (Statistical Coupling analysis) - Lockless and Ranganathan, 1999.
    As on Halperin et al, 2006.
    * ELSC (Explicit Likelihood of Subset Covariation) - Dekker et al, 2004.
    """
    
    def __init__(self, file1, file2, id1, id2, chain1, chain2, 
                 alignment, coevolution):
        self.file1 = str(file1)
        self.file2 = str(file2)
        self.id1 = str(id1)
        self.id2 = str(id2)
        self.chain1 = str(chain1)
        self.chain2 = str(chain2)
        self.alignment = alignment
        self.coevolution = coevolution
        
    def __call__(self, file1, file2, id1, id2, chain1, chain2, 
                 alignment, coevolution):
        self.file1 = str(file1)
        self.file2 = str(file2)
        self.id1 = str(id1)
        self.id2 = str(id2)
        self.chain1 = str(chain1)
        self.chain2 = str(chain2)
        self.alignment = alignment
        self.coevolution = coevolution

    def coevolAnalysis(self, id1, id2, chain1, chain2, alignment, coevolution):
        "Returns a matrix of coevolution scores"
        
        seq = class_sequence(self.file1, self.file2, self.id1, self.id2, 
                       self.chain1, self.chain2)
        aln = class_alignment(self.id1, self.id2, self.alignment)
        
        alignment1 = aln.cutAlignment(id1, alignment)
        alignment2 = aln.cutAlignment(id2, alignment)
        
        protein1 = []
        protein2 = []
        try:
            protein1 = seq.matchResiduePosition(id1, chain1)
            protein2 = seq.matchResiduePosition(id2, chain2)
        except:
            pass

        info = dict()
        
        if coevolution == "mi":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)
            pD1 = probabilityDict(columns1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            pD2 = probabilityDict(columns2)
         
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    info[(i,j)] = mutualInformation(i, j, columns1, columns2, pD1, pD2)
        
        elif coevolution == "mie":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)
            pD1 = probabilityDict(columns1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            pD2 = probabilityDict(columns2)
         
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    info[(i,j)] = miEntropy(i, j, columns1, columns2, pD1, pD2)
                    
        elif coevolution == "rcwmi":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)
            pD1 = probabilityDict(columns1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            pD2 = probabilityDict(columns2)
         
            i_all = dict()
            all_j = dict()
            for i in range(len(columns1)):
                v_i = 0
                for j in range(len(columns2)):
                    v_i += mutualInformation(i, j, columns1, columns2, 
                                             pD1, pD2)
                    i_all[i]= v_i

            for j in range(len(columns2)):
                v_j = 0
                for i in range(len(columns1)):
                    v_j += mutualInformation(i, j, columns1, columns2, 
                                             pD1, pD2)
                    all_j[j]= v_j
            
            column = columns1[0]
            n = len(column)
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    mi = mutualInformation(i, j, columns1, columns2, 
                                           pD1, pD2)    
                    info[(i,j)] = rowColumnWeighed(mi, 
                                                   i_all[i], all_j[j], n)
        
        elif coevolution == "cpvnmie":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)
            pD1 = probabilityDict(columns1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            pD2 = probabilityDict(columns2)
            
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    res1 = str(alignment1[0][i])
                    res2 = str(alignment2[0][j])
                    mie = miEntropy(i, j, columns1, columns2, pD1, pD2)
                    info[(i,j)] = contactPreferenceMI(mie, res1, res2)
                    
        elif coevolution == "cpvn":
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    average = []
                    for a,b in zip(columns1[i],columns2[j]):
                        if a in aa and b in aa:
                            average.append(float(matchScore(res1, res2, "CPVN")))
                    info[(i,j)] = mean(average)

        elif coevolution == "clm":
            for i in range(len(alignment1[0])):
                for j in range(len(alignment2[0])):
                    average = []
                    for a,b in zip(columns1[i],columns2[j]):
                        if a in aa and b in aa:
                            average.append(float(matchScore(res1, res2, "CLM")))
                    info[(i,j)] = mean(average)
                    
        elif coevolution == "vol":
            for i in range(len(alignment1[0])):
                for j in range(len(alignment2[0])):
                    average = []
                    for a,b in zip(columns1[i],columns2[j]):
                        if a in aa and b in aa:
                            average.append(float(matchScore(res1, res2, "VOL")))
                    info[(i,j)] = mean(average)
                    
        elif coevolution == "omes":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            omes = dict()
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    omes[(i,j)] = covarianceOMES(columns1[i],columns2[j])
            max_pos = []
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    max_pos.append(omes[(i,j)])
            max_val = max(max_pos)
                    
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if omes[(i,j)] != 0.0:
                        info[(i,j)] = omes[(i,j)] * 1.0 / max_val
                    else:
                        info[(i,j)] = 0.0
                    
        elif coevolution == "pearson":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            score_matrix = mapMatrix("MCLACHLAN")
            N = len(columns1[0])
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    d_matrix1 = twoDimensionalMatrix(columns1[i], score_matrix)
                    d_matrix2 = twoDimensionalMatrix(columns2[j], score_matrix)
                    info[(i,j)] = pearsonsCorrelation(d_matrix1, d_matrix2, N)
                    
        elif coevolution == "spearman":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            score_matrix = mapMatrix("MCLACHLAN")
            spearman = dict()
            N = len(columns1[0])
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    d_matrix1 = twoDimensionalMatrix(columns1[i], score_matrix)
                    d_matrix2 = twoDimensionalMatrix(columns2[j], score_matrix)
                    spearman[(i,j)] = spearmansCorrelation(d_matrix1, d_matrix2, N)
            
            max_pos = []
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    max_pos.append(spearman[(i,j)])
            max_val = max(max_pos)
                    
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if spearman[(i,j)] != 0.0:
                        info[(i,j)] = spearman[(i,j)] * 1.0 / max_val
                    else:
                        info[(i,j)] = 0.0
                    
        elif coevolution == "mcbasc":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            score_matrix = mapMatrix("MCLACHLAN")
            N = len(columns1[0])
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    d_matrix1 = twoDimensionalMatrix(columns1[i], score_matrix)
                    d_matrix2 = twoDimensionalMatrix(columns2[j], score_matrix)
                    info[(i,j)] = mcbascCorrelation(d_matrix1,d_matrix2, N)
                     
        
        elif coevolution == "quartets":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            quartets = dict()
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    quartets[(i,j)] = quartetsCorrelation(columns1[i],columns2[j])
            
            max_pos = []
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    max_pos.append(quartets[(i,j)])
            max_val = max(max_pos)
                    
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if quartets[(i,j)] != 0.0:
                        info[(i,j)] = quartets[(i,j)] * 1.0 / max_val
                    else:
                        info[(i,j)] = 0.0
                        
        elif coevolution == "sca":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
            
            sca = dict()   
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    sca[(i,j)] = perturbationSCA(columns1[i],columns2[j],\
                                                  j,columns2)
            max_pos = []
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    max_pos.append(sca[(i,j)])
            max_val = max(max_pos)
                    
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if sca[(i,j)] != 0.0:
                        info[(i,j)] = sca[(i,j)] * 1.0 / max_val
                    else:
                        info[(i,j)] = 0.0
                    
        elif coevolution == "elsc":
            alignment1 = [e for e in alignment1]
            columns1 = transpose(alignment1)

            alignment2 = [e for e in alignment2]
            columns2 = transpose(alignment2)
             
            elsc = dict()  
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    elsc[(i,j)] = perturbationELSC(columns1[i],columns2[j],\
                                                   j,columns2)       
            max_pos = []
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    max_pos.append(elsc[(i,j)])
            max_val = max(max_pos)
                    
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if elsc[(i,j)] != 0.0:
                        info[(i,j)] = elsc[(i,j)] * 1.0 / max_val
                    else:
                        info[(i,j)] = 0.0               
        else: pass
        
        
        output = "./Results/" + alignment + "_" + coevolution + ".txt"
        results = open(output, "w")
        for i, j in sorted(info.keys()):
            if protein1 != [] and protein2 != []:
                print >> results, protein1[i], protein2[j], \
                round((info[(i, j)]), 4)
            elif protein1 != [] and protein2 == []:
                print >> results, protein1[i], protein1[j], \
                round((info[(i, j)]), 4)
            else:
                print >> results, str(i+1), str(j+1), \
                round((info[(i, j)]), 4)
        results.close()
    
    def bestInfo(self, id1, id2, alignment, coevolution):
        "Points out the best coevolution scores"
        
        seq = class_sequence(self.file1, self.file2, self.id1, self.id2, 
                       self.chain1, self.chain2)
        
        histogram = results_histogram
        heatmap = results_heatmap
        best_info = best_results
        
        surface1 = []
        surface2 = []
        interface = []
        try:
            surface1 = seq.parseSurfacePDB(id1)
            surface2 = seq.parseSurfacePDB(id2)
        except:
            pass
        
        try:
            interface = seq.parseInterfacePDB(id1)
        except:
            pass
        
        input = "./Results/" + alignment + "_" + coevolution + ".txt"
        output = "./Results/" + alignment + "_" + coevolution + "_best.txt"
        bestResults(input, output, best_info, surface1, surface2, interface)
        
        if histogram == True:
            input = "./Results/" + alignment + "_" + coevolution + ".txt"
            output = "./Results/" + alignment + "_" + coevolution + "_hg.png"
            drawHistogram(input, output)
            
        if heatmap == True:
            input = "./Results/" + alignment + "_" + coevolution + ".txt"
            output = "./Results/" + alignment + "_" + coevolution + "_hm.png"
            drawHeatmap(id1, id2, input, output)
        
        
    def structureSingle(self, id1, id2, chain1, chain2, alignment, coevolution):
        "Structure based results for proteins with single chain"
        
        structure = results_structure
        best_info = best_results
        
        input = "./Results/" + alignment + "_" + coevolution + "_best.txt"
        input_results = open(input, "r")
        results = input_results.readlines()
        input_results.close()
            
        positions1 = []
        positions2 = []
        for line in results:
            l = line.rstrip("\n")
            l = l.split()
            res1 = int(l[0])
            res2 = int(l[1])
            positions1.append(res1)
            positions2.append(res2)
            
        if structure == "pymol":
            output1 = "./Data/" + id1 + ".pml"
            out_struct1 = open(output1, "w")
            print >> out_struct1, "load %s" %(id1 + ".pdb")
            print >> out_struct1, "hide lines"
            print >> out_struct1, "hide nonbonded"
            print >> out_struct1, "bg_color black"
            print >> out_struct1, "color grey20"
            print >> out_struct1, "show cartoon"
            print >> out_struct1, "select hitmol, chain %s" %(chain1.lower())
            print >> out_struct1, "color red, (hitmol and resid *)"
            for pos in positions1:
                if len(positions1) <= 20:
                    print >> out_struct1, "color yellow, (hitmol and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct1, "show spheres, (hitmol and resid %s)" \
                    %(str(pos +1))     
                else:
                    print >> out_struct1, "color yellow, (hitmol and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct1, "show sticks, (hitmol and resid %s)" \
                    %(str(pos +1))
            out_struct1.close()
            
            output2 = "./Data/" + id2 + ".pml"
            out_struct2 = open(output2, "w")
            print >> out_struct2, "load %s" %(id2 + ".pdb")
            print >> out_struct2, "hide lines"
            print >> out_struct2, "hide nonbonded"
            print >> out_struct2, "bg_color black"
            print >> out_struct2, "color grey20"
            print >> out_struct2, "show cartoon"
            print >> out_struct2, "select hitmol, chain %s" %(chain2.lower())
            print >> out_struct2, "color blue, (hitmol and resid *)"
            for pos in positions2:
                if best_info <= 20:
                    print >> out_struct2, "color green, (hitmol and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct2, "show spheres, (hitmol and resid %s)" \
                    %(str(pos +1))     
                else:
                    print >> out_struct2, "color green, (hitmol and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct2, "show sticks, (hitmol and resid %s)" \
                    %(str(pos +1))
            out_struct2.close()
        else: pass

        
    def structurePair(self, id1, id2, chain1, chain2, alignment, coevolution):
        "Structure based results for a protein with two chains"
        
        structure = results_structure
        best_info = best_results
        
        input = "./Results/" + alignment + "_" + coevolution + "_best.txt"
        input_results = open(input, "r")
        results = input_results.readlines()
        input_results.close()
            
        positions1 = []
        positions2 = []
        for line in results:
            l = line.rstrip("\n")
            l = l.split()
            res1 = int(l[0])
            res2 = int(l[1])
            positions1.append(res1)
            positions2.append(res2)
            
        if structure == "pymol":
            output = "./Data/" + id1 + ".pml"
            
            out_struct = open(output, "w")
            print >> out_struct, "load %s" %(id1 + ".pdb")
            print >> out_struct, "hide lines"
            print >> out_struct, "hide nonbonded"
            print >> out_struct, "bg_color black"
            print >> out_struct, "color grey20"
            print >> out_struct, "show cartoon"
            print >> out_struct, "select hitmol1, chain %s" %(chain1.lower())
            print >> out_struct, "select hitmol2, chain %s" %(chain2.lower())
            print >> out_struct, "color red (hitmol1)"
            print >> out_struct, "color blue (hitmol2)" + "\n"
            for pos in positions1:
                if best_info <= 20:
                    print >> out_struct, "color yellow, (hitmol1 and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct, "show spheres, (hitmol1 and resid %s)" \
                    %(str(pos +1))     
                else:
                    print >> out_struct, "color yellow, (hitmol1 and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct, "show sticks, (hitmol1 and resid %s)" \
                    %(str(pos +1))
                    
            for pos in positions2:
                if best_info <= 20:
                    print >> out_struct, "color green, (hitmol2 and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct, "show spheres, (hitmol2 and resid %s)" \
                    %(str(pos +1))     
                else:
                    print >> out_struct, "color green, (hitmol2 and resid %s)" \
                    %(str(pos +1))
                    print >> out_struct, "show sticks, (hitmol2 and resid %s)" \
                    %(str(pos +1))    
            out_struct.close()
        else: 
            pass      

        
def matchScore(alpha, beta, score_matrix):
    "Matches scores from a matrix"
        
    alphabet = {}    
    alphabet["I"] = 0
    alphabet["V"] = 1
    alphabet["L"] = 2
    alphabet["F"] = 3
    alphabet["C"] = 4
    alphabet["M"] = 5
    alphabet["A"] = 6
    alphabet["G"] = 7
    alphabet["T"] = 8
    alphabet["S"] = 9
    alphabet["W"] = 10
    alphabet["Y"] = 11
    alphabet["P"] = 12
    alphabet["H"] = 13
    alphabet["E"] = 14
    alphabet["Q"] = 15
    alphabet["D"] = 16
    alphabet["N"] = 17
    alphabet["K"] = 18
    alphabet["R"] = 19
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return score_matrix[lut_x][lut_y]

def matchScore2(alpha, beta, score_matrix):
    "Matches scores from a matrix - different residue order"
    
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
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return score_matrix[lut_x][lut_y]
    
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
   
def twoDimensionalMatrix(column, score_matrix):
    "For each column in the alignment constructs a two-dimensional matrix"
    
    two_d = []
    for i in range(len(column)):
        for j in range(len(column)):
            if i!=j:
                res1 = column[i]
                res2 = column[j]
                if res1 in aa and res2 in aa:
                    s = float(matchScore2(res1, res2, score_matrix))
                    two_d.append(s)
                else:
                    s = 0.0
                    two_d.append(s)
                
    return two_d
        
def log21(n):  
    return log(n) * 1.0 / log(21)

def ln(n): 
    return log(n) * 1.0 / log(e)

def transpose(L):
    R = range(len(L[0]))
    rL = list()
    for i in R:
        rL.append(''.join([item[i] for item in L]))
    return rL
    

def probabilityDict(columns):
    "Caches character probabilities for each column"
    
    n = len(columns[0])
    pD = list()
    for col in columns:
        aa = list(set(col))
        values = [col.count(k) * 1.0 / n for k in aa]
        pD.append(dict(zip(aa, values)))
    return pD


def mutualInformation(i, j, cols1, cols2, pD1, pD2):
    """
    Mutual informaton for protein coevolution as by
    Gloor et al, 2005. MI(X,Y) = H(X) + H(Y) - H(X,Y)
    MI(X,Y) = SUMSUM P(x,y).log20(P(x,y)/P(x).P(y))
    Treates gaps as signal.
    """
    
    col1, col2 = cols1[i], cols2[j]
    n = len(col1)
    assert n == len(col2)
    mi = 0
    pairs = [col1[k] + col2[k] for k in range(n)]
    pL = sorted(list(set(pairs)))
    for p in pL: 
        pXY = pairs.count(p) * 1.0 / n
        pX = pD1[i][p[0]]
        pY = pD2[j][p[1]]
        inside = (pXY * 1.0) / (pX * pY)
        outside = pXY * log21(inside)
        mi += outside
    return mi

def miEntropy(i, j, cols1, cols2, pD1, pD2):
    """
    Mutual informaton by pair entropy - Martin et al, 2005.
    MI(X,Y) = (H(X) + H(Y) - H(X,Y)) / H(X,Y)
    MI(X,Y) = (SUMSUM P(x,y).log20(P(x,y)/P(x).P(y))) / 
               -(SUMSUM P(x,y).log20(P(x,y)))  
    """

    col1, col2 = cols1[i], cols2[j]
    assert len(col1) == len(col2)
    n = len(col1)
    mi = 0
    entropy = 0
    pairs = [col1[k] + col2[k] for k in range(n)]
    pL = sorted(list(set(pairs)))
    for p in pL: 
        pXY = pairs.count(p) * 1.0 / n
        pX = pD1[i][p[0]]
        pY = pD2[j][p[1]]
        inside = (pXY * 1.0) / (pX * pY)
        outside = pXY * log21(inside)
        mi += outside
    for p in pL: 
        pXY = pairs.count(p) * 1.0 / n
        inside = pXY
        outside = pXY * log21(inside)
        entropy += outside
    entropy = -entropy
    if entropy == 0.0:
        mi_entropy = 0.0
    else: mi_entropy = mi / entropy
    return mi_entropy
  
def rowColumnWeighed(mi, i_all, all_j, n):
    """
    Row and Column weighed Mutual Information - Gouveia-
    Oliveira et al, 2007. 
    RCW(X,Y) = MI(X,Y) / 
            (((MI(X,all) + MI(all,Y) - 2MI(X,Y))/(n-1))
    """
    
    bottom = (i_all + all_j - 2.0 * mi) / (n - 1)
    if bottom == 0.0:
        rcwmi = 0.0
    else: rcwmi = mi / bottom
    
    return rcwmi

def contactPreferenceMI(mie, res1, res2):
    """
    Contact preferences, volume normalized MIE - F. Madeira et al, 2012
    CPVN MI/E(X,Y) = k/2 * CPVN(X,Y) + MI(X,Y)/2H(X,Y);
    """
    
    k = 0.1 
    cpvn = float(matchScore(res1, res2, "CPVN")) 
    cpvnmie = (k / 2.0 * cpvn) + (mie / 2.0)
    
    return cpvnmie

def covarianceOMES(column1,column2):
    """
    Normalized Covariance analysis; OMES - Observed Minus Expected Squared
    derived from the covariance method of Kass and Horovitz, 2002
    """
 
    assert len(column1) == len(column2)
    
    L = []
    Nvalid = []
    Cxi = []
    Cyj = []
    for i,j in zip(column1,column2):
        if i in aa and j in aa:
            value = [i,j]
            Nvalid.append(value)
            Cxi.append(i)
            Cyj.append(j)
            if value not in L:
                L.append(value)

    len_Nvalid = len(Nvalid)
    omes = 0.0
    for value in L:
        Nobs = Nvalid.count(value)
        i = value[0]
        j = value[1]
        Ci = Cxi.count(i)
        Cj = Cyj.count(j)
        Nex = Ci * Cj / len_Nvalid    
        top = (Nobs - Nex)**2
        omes += top * 1.0 / len_Nvalid
    
    return omes

def pearsonsCorrelation(d_matrix1,d_matrix2, N):
    """
    Pearson's Correlation (Gobel method) - Gobel et al, 1994.
    """
    
    assert len(d_matrix1) == len(d_matrix2)
    
    no_match = 0.0
    for k,l in zip(d_matrix1,d_matrix2):
        if k!=l:
            no_match += 1.0
    length = len(d_matrix1)
    Wkl = no_match * 1.0 / length
    
    sigma_i = std(d_matrix1)
    Si = []
    av_Si = mean(d_matrix1)
    for i in (d_matrix1):
        Si.append(i - av_Si)
    
    sigma_j = std(d_matrix2)
    Sj = []
    av_Sj = mean(d_matrix1)
    for j in (d_matrix2):
        Sj.append(j - av_Sj)
    
    top = 0.0
    for i,j in zip(Si,Sj):
        top += float(i * j * Wkl)

    bottom = sigma_i * sigma_j
    if bottom == 0.0:
        pearson = 0.0
    else:
        pearson = (1.0 / N**2)*(top/bottom)
    
    return pearson

def spearmansCorrelation(d_matrix1,d_matrix2, N):
    """
    Spearman's rank Correlation - Pazos et al, 1997. 
    """
    
    assert len(d_matrix1) == len(d_matrix2)
    
    rank_matrix1 = []
    rank_matrix2 = []
    rank_temp1 = []
    rank_temp2 = []
    for k,l in zip(d_matrix1,d_matrix2):
        if k not in rank_temp1:
            rank_temp1.append(k)
            cnt = d_matrix1.count(k)
            rank = cnt * 1.0 / len(d_matrix1)
            rank_matrix1.append(rank)
        if l not in rank_temp2:
            rank_temp2.append(l)
            cnt = d_matrix2.count(l)
            rank = cnt * 1.0 / len(d_matrix2)
            rank_matrix2.append(rank)
    
    no_match = 0.0
    for k,l in zip(d_matrix1,d_matrix2):
        if k!=l:
            no_match += 1.0
    length = len(d_matrix1)
    Wkl = no_match * 1.0 / length
    
    sigma_i = std(d_matrix1)
    Si = []
    av_Si = mean(d_matrix1)
    for i in (rank_matrix1):
        Si.append(i - av_Si)
    
    sigma_j = std(d_matrix2)
    Sj = []
    av_Sj = mean(d_matrix1)
    for j in (rank_matrix2):
        Sj.append(j - av_Sj)
    
    top = 0.0
    for i,j in zip(Si,Sj):
        top += float(i * j * Wkl)

    bottom = sigma_i * sigma_j
    if bottom == 0.0:
        spearman = 0.0
    else:
        spearman = (1.0 / N**2)*(top/bottom)
    
    return spearman

def mcbascCorrelation(d_matrix1,d_matrix2, N):
    """
    McBASC - McLachlan Based Substitution Correlation.
    Fodor and Aldrich, 2004.
    """
    
    assert len(d_matrix1) == len(d_matrix2)
    
    sigma_i = std(d_matrix1)
    Si = []
    av_Si = mean(d_matrix1)
    for i in (d_matrix1):
        Si.append(i - av_Si)
    
    sigma_j = std(d_matrix2)
    Sj = []
    av_Sj = mean(d_matrix1)
    for j in (d_matrix2):
        Sj.append(j - av_Sj)
    
    top = 0.0
    for i,j in zip(Si,Sj):
        top += float(i * j)

    bottom = sigma_i * sigma_j
    if bottom == 0.0:
        mcbasc = 0.0
    else:
        mcbasc = abs((1.0 / N**2)*(top/bottom))
    
    return mcbasc


def quartetsCorrelation(column1,column2):
    """
    Normalized Quartets correlation method by Galitsky, 2002.
    """
 
    assert len(column1) == len(column2)
    
    quartets = 0.0
    x = column1
    y = column2
    pairs = []
    for i,j in zip(x,y):
        value = [i,j]
        pairs.append(value)
        
    for i,j in zip(x,y):
        if i in aa and j in aa:
            Pix = x.count(i) * 1.0 / len(x)
            Piy = y.count(i) * 1.0 / len(y)
            Pjx = x.count(j) * 1.0 / len(x)
            Pjy = y.count(j) * 1.0 / len(y)
            val = [i,j]
            Dmin = pairs.count(val)
            Dif = 1.0 * (len(pairs) - Dmin)
            if Dif != 0.0:
                DQmin = Dmin * 1.0 / Dif
            else:
                DQmin = 0.0

            try :
                if ((Pix*Pjy > Piy*Pjx) and ((Pix > Dmin) or (Pjy > Dmin)) or\
                    (Pix*Pjy < Piy*Pjx) and ((Piy > Dmin) or (Pjx > Dmin)))\
                    and\
                   (((Pix*Pjy) * 1.0 / (Piy*Pjx) > DQmin) or\
                    ((Piy*Pjx) * 1.0 / (Pix*Pjy) > DQmin)):
                    quartets += 1.0
            except:
                quartets += 0
    return quartets

def perturbationSCA(column1, column2, j, columns2):
    """
    Normalized SCA - Statistical Coupling analysis, Lockless and 
    Ranganathan, 1999. As on Halperin et al, 2006.
    """
    
    assert len(column1) == len(column2)
    
    new_columns2 = subAlignment(column2, columns2)
    x = column1
    y = new_columns2[j]
    
    inside = 0.0
    for i in x:
        if i in aa:
            Pix = x.count(i) * 1.0 / len(x)
            Pixj = y.count(i) * 1.0 / len(y)
            if Pixj != 0.0:
                inside += (ln(Pixj) - Pix)**2
            
    sca = sqrt(inside)
    return sca

def perturbationELSC(column1, column2, j, columns2):
    """
    Normalized ELSC - Explicit Likelihood of Subset Covariation, 
    Dekker et al, 2004.
    """
    
    assert len(column1) == len(column2)
    
    new_columns2 = subAlignment2(column1, column2, columns2)
    x = column1
    y1= column2
    y2 = new_columns2[j]
    
    
    comb_x = []
    comb_all = []
    for i in x:
        if i in aa:
            Nxj = y1.count(i)
            nxj = y2.count(i)
            Nall = len(y1)
            nall = len(y2)
            mxj = int(round((Nxj * 1.0 / Nall) * nall))
            top = long(factorial(Nxj))
            bot1 = factorial(nxj)* factorial(Nxj - nxj)
            bot2 = factorial(mxj)* factorial(Nxj - mxj)
            comb_x.append(top / bot1)
            comb_all.append(top / bot2)          
    
    product = 1.0
    for k,l in zip(comb_x, comb_all):    
        product *= (k * 1.0 / l) 
        
    if product != 0.0:
        elsc = -ln(product)
    else: 
        elsc = 0.0
    
    return elsc

def subAlignment (column, columns):
    "Creates a sub_alignment based on the most frequent AA in column"
    
    pD = []
    y = column
    for j in range(len(y)):
        if y[j] in aa:
            freq = y.count(y[j])
            freq_aa = y[j]
            value = [freq_aa, freq]
            pD.append(value)
    
    sort = sorted(pD, key=lambda pD: pD[1])
    aa_j = sort[0][0]
    
    col_positions = []
    pos = -1
    for j in y:
        pos += 1
        if j == aa_j:
            col_positions.append(pos)
    
    sub_align = []
    for col in columns:
        sub_col = []
        for pos in col_positions:
            sub_col.append(col[pos])
        sub_align.append(sub_col)
    return sub_align

def subAlignment2 (column1, column2, columns):
    "Creates a sub_alignment based on AA identity of column1"
    
    x = column1
    y = column2
    
    list_i = []
    for i in x:
        if i in aa:
            if i not in list_i:
                list_i.append(i)
    
    col_positions = []
    pos = -1
    for j in y:
        pos += 1
        if j in list_i:
            col_positions.append(pos)
    
    sub_align = []
    for col in columns:
        sub_col = []
        for pos in col_positions:
            sub_col.append(col[pos])
        sub_align.append(sub_col)
    return sub_align
  
def bestResults(input, output, best_info, surface1, surface2, interface):
    "Creates a new list of best coevolution scores"
    
    input_results = open(input, "r")
    results = input_results.readlines()
    input_results.close()
    
    all = []
    for line in results:
        if line == "\n": pass
        else: 
            l = line.rstrip("\n")
            l = l.split()
            res1 = int(l[0])
            res2 = int(l[1])
            mi = float(l[2])
            if res1 in surface1 and res2 in surface2:
                value = [res1, res2, mi]
                all.append(value)
            elif res1 in surface1 and res2 in surface1:
                value = [res1, res2, mi]
                all.append(value)
            else:
                value = [res1, res2, mi]
                all.append(value)
    
    a = all
    sort = sorted(a, key=lambda a: a[2])
    length = len(sort)
    position = length - best_info
    threshold = sort[position]
    
    out_best = open(output, "w")
    for line in all:
        res1 = line[0]
        res2 = line[1]
        mi = float(line[2])
        value = [res1, res2]
        if mi >= threshold[2]:
            if value in interface:
                print >> out_best, res1, res2, mi, "Interface contact"
            else:
                print >> out_best, res1, res2, mi
    out_best.close()

def drawHistogram(input, output):
    "Creates a histogram of coevolution scores"  
        
    data = []
    info = []
    input_results = open(input, "r")
    results = input_results.readlines()
    input_results.close()
        
    for line in results:
        l = line.rstrip("\n")
        l = l.split()
        res1 = int(l[0])
        res2 = int(l[1])
        inf = float(l[2])
        value = [res1, res2, inf]
        data.append(value)
        info.append(inf)

    maxi = max(info)
    L = [t[2] for t in data]
    X = maxi
    pyplot.hist(L, bins=X * 50)
    ax = pyplot.axes()
    ax.set_xlabel('Score')
    ax.set_ylabel('Frequency')
    ax.set_xlim(0, X)
    pyplot.savefig(output)
                
def drawHeatmap(id1, id2, input, output): 
    "Creates a heatmap of coevolution scores"
               
    input_results = open(input, "r")
    results = input_results.readlines()
    input_results.close()
    
    data = []  
    residue1 = []
    residue2 = [] 
    for line in results:
        l = line.rstrip("\n")
        l = l.split()
        res1 = int(l[0])
        res2 = int(l[1])
        inf = float(l[2])
        value = [res1, res2, inf]
        data.append(value)
        residue2.append(res1)
        if res1 not in residue1:
            residue1.append(res1)
       
    startX = int(data[0][0])
    startY = int(data[0][1])
    length = len(data)
    endX = int(data[length -1][0])
    endY = int(data[length -1][1])
    
    lenX = len(residue1)
    lenY = residue2.count(startX)    
        
    heatmap = zeros((lenX+1, lenY+1))
    for i in range(0, len(data)-1):
        X = int(data[i][0])
        Y = int(data[i][1])
        XY = float(data[i][2])
        heatmap[X][Y] = XY
            
    pyplot.figure()
    pyplot.pcolormesh(heatmap)
    pyplot.colorbar() 
    pyplot.axes().set_xlabel(id1)
    pyplot.axes().set_ylabel(id2)
    pyplot.axes().set_xlim(startX, endX)
    pyplot.axes().set_ylim(startY, endY)
    pyplot.savefig(output)


