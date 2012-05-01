###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

"""Parameters Loader"""

from ConfigParser import SafeConfigParser

surface_threshold = 7                   # 7% [0, max(surface)[ 
psiblast_evalue = 0.00001               # [0.00000001:10]
psiblast_identity = 30                  # [0:100] (%)
psiblast_coverage = 60                  # [0:100] (%)
psiblast_threading = False              # Number of cores/servers or False
pairwise_distance = "clustalw"          # "clustalw", "pdistance", "Kimura"
                                        # "jukescantor" or "alignscore"
alignscore_matrix = "BLOSUM62"          # "BLOSUM62" or "PAM250"
theilsen_cutoff = 0.7                   # [0.25:1.0(all sequences)]
clustalw_gap_opening = 10               # [0:100]
clustalw_gap_extension = 0.2            # [0:10]
clustalw_distance_matrix = "GONNET"     # "GONNET", "BLOSUM" or "PAM"
muscle_max_iteration = 16               # [2:16]
mafft_configuration = "linsi"           # "fftnsi" or "linsi"
mafft_threading = False                 # Number of cores/servers or False
alphabet_reduction = False              # False or "charge", "charge_his", "polarity" 
                                        # or "hydropathy"
alignment_score = False                 # "sumofpairs" or False
best_results = 20                       # [1:max(scores)]
results_histogram = True                # True or False
results_heatmap = True                  # True or False
results_structure = "pymol"             # "pymol" or False
results_sifts = False                   # True or False

def SaveParameters(filename):
    "Saves default parameters"
    
    parser = SafeConfigParser()
    parser.add_section('Global')
    parser.add_section('Psiblast')
    parser.add_section('Clustalw')
    parser.add_section('Muscle')
    parser.add_section('Mafft')
    parser.add_section('Results')
    parser.set('Global', 'SurfaceThreshold', float(surface_threshold))
    parser.set('Psiblast', 'Evalue', float(psiblast_evalue))
    parser.set('Psiblast', 'Identity', int(psiblast_identity))
    parser.set('Psiblast', 'Coverage', int(psiblast_coverage))
    parser.set('Psiblast', 'Threading', str(psiblast_threading))
    parser.set('Global', 'PairwiseDistance', str(pairwise_distance))
    parser.set('Clustalw', 'GapOpening', float(clustalw_gap_opening))
    parser.set('Clustalw', 'GapExtension', float(clustalw_gap_extension))
    parser.set('Clustalw', 'Matrix', str(clustalw_distance_matrix))
    parser.set('Global', 'Matrix', str(alignscore_matrix))
    parser.set('Global', 'TheilSenCutoff', float(theilsen_cutoff))
    parser.set('Muscle', 'MaxIteration', int(muscle_max_iteration))
    parser.set('Mafft', 'Configuration', str(mafft_configuration))
    parser.set('Mafft', 'Threading', str(mafft_threading))
    parser.set('Global', 'AlphabetReduction', str(alphabet_reduction))
    parser.set('Global', 'AlignmentScore', str(alignment_score))
    parser.set('Results', 'Best', int(best_results))
    parser.set('Results', 'Histogram', str(results_histogram))
    parser.set('Results', 'Heatmap', str(results_heatmap))
    parser.set('Results', 'Structure', str(results_structure))
    parser.set('Results', 'Sifts', str(results_sifts))
    fil = open(filename, 'w')
    parser.write(fil)
    fil.close()

def LoadParameters(filename, option):
    "Loads and tests input parameters"
    
    parser = SafeConfigParser()
    try:
        parser.read(filename)
        if option == "surface_threshold":
            surface_threshold = parser.getfloat('Global', 'SurfaceThreshold') 
            return surface_threshold
        elif option == "psiblast_evalue":
            psiblast_evalue = parser.getfloat('Psiblast', 'Evalue')
            return psiblast_evalue
        elif option == "psiblast_identity":
            psiblast_identity = parser.getint('Psiblast', 'Identity')
            return psiblast_identity
        elif option == "psiblast_coverage":
            psiblast_coverage = parser.getint('Psiblast', 'Coverage')
            return psiblast_coverage
        elif option == "psiblast_threading":
            psiblast_threading = parser.get('Psiblast', 'Threading')
            return psiblast_threading
        elif option == "pairwise_distance":
            pairwise_distance = parser.get('Global', 'PairwiseDistance')
            return pairwise_distance
        elif option == "clustalw_gap_opening":
            clustalw_gap_opening = parser.getfloat('Clustalw', 'GapOpening')
            return clustalw_gap_opening
        elif option == "clustalw_gap_extension":
            clustalw_gap_extension = parser.getfloat('Clustalw', 'GapExtension')
            return clustalw_gap_extension
        elif option == "clustalw_distance_matrix":
            clustalw_distance_matrix = parser.get('Clustalw', 'Matrix')
            return clustalw_distance_matrix
        elif option == "alignscore_matrix":
            alignscore_matrix = parser.get('Global', 'Matrix')
            return alignscore_matrix
        elif option == "theilsen_cutoff":
            theilsen_cutoff = parser.getfloat('Global', 'TheilSenCutoff')
            return theilsen_cutoff
        elif option == "muscle_max_iteration":
            muscle_max_iteration = parser.getint('Muscle', 'MaxIteration')
            return muscle_max_iteration
        elif option == "mafft_configuration":
            mafft_configuration = parser.get('Mafft', 'Configuration')
            return mafft_configuration
        elif option == "mafft_threading":
            mafft_threading = parser.get('Mafft', 'Threading')
            return mafft_threading
        elif option == "alphabet_reduction":
            alphabet_reduction = parser.get('Global', 'AlphabetReduction')
            return alphabet_reduction
        elif option == "alignment_score":
            alignment_score = parser.get('Global', 'AlignmentScore')
            return alignment_score
        elif option == "best_results":
            best_results = parser.getint('Results', 'Best')
            return best_results
        elif option == "results_histogram":
            results_histogram = parser.getboolean('Results', 'Histogram')
            return results_histogram
        elif option == "results_heatmap":
            results_heatmap = parser.getboolean('Results', 'Heatmap')
            return results_heatmap
        elif option == "results_structure":
            results_structure = parser.get('Results', 'Structure')
            return results_structure
        elif option == "results_sifts":
            results_sifts = parser.getboolean('Results', 'Sifts')
            return results_sifts
        elif option == "test":
            parser.getint('Results', 'Best')
            print "Parameters... OK"
            return
        else:
            raise StandardError, "ERROR: Invalid option"
    except:
        raise StandardError, "ERROR: Invalid Parameters File"

