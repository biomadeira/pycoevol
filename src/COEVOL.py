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
from math import log
from numpy import mean, zeros
from matplotlib import pyplot


class coevolution:
    """
    Main code for coevolution analysis.
    
    Matrix-based Methods:
    Residue Contact Preferences, Volume Normalized - Glaser et al, 2001.
    Contact likelihood matrix - Singer et al, 2002.
    Residue volume, normalized  - Esque et al, 2010
    
    Mutual Information based methods:
    Mutual Information - Martin el al, 2005
    MI by pair entropy - Martin el al, 2005
    Row and column weighed MI - Gouveia-Oliveira et al, 2007
    Contact preferences, volume normalized MIE - unpublished
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
         
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    mi = mutualInformation(i, j, columns1, columns2, 
                                           pD1, pD2)
                    
                    i_all = 0
                    for Js in range(len(columns2)):
                        i_all += mutualInformation(i, Js, 
                                                   columns1, columns2, pD1, pD2)
                    
                    all_j = 0
                    for Is in range(len(columns1)):
                        all_j += mutualInformation(Is, j, 
                                                   columns1, columns2, pD1, pD2)
                    
                    column = columns1[0]
                    n = len(column)        
                    info[(i,j)] = rowColumnWeighed(mi, i_all, all_j, n)
        
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
        else: pass
        
        
        output = "./Data/" + alignment + "_" + coevolution + ".txt"
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
        
        input = "./Data/" + alignment + "_" + coevolution + ".txt"
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

        
def matchScore(alpha, beta, matrix):
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
   
    
def log20(n):  
    return log(n) * 1.0 / log(20)

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
    Martin et al, 2005. MI(X,Y) = H(X) + H(Y) - H(X,Y)
    MI(X,Y) = SUMSUM P(x,y).log20(P(x,y)/P(x).P(y))
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
        outside = pXY * log20(inside)
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
        outside = pXY * log20(inside)
        mi += outside
    for p in pL: 
        pXY = pairs.count(p) * 1.0 / n
        inside = pXY
        outside = pXY * log20(inside)
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


