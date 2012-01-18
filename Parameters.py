###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

"Editable parameters"

surface_threshold = 7                   # 7% (Angstrom^2), [0, max(surface)] 
interface_threshold = 1                 # 1 Angstrom^2, [0, max(surface)]
psiblast_evalue = 10                    # [0.0000001:10]
psiblast_identity = 0                   # [0:100] (%)
psiblast_coverage = 0                   # [0:100] (%)
pairwise_trim = True                    # True or False
pairwise_distance = "jukescantor"       # "pdistance", "jukescantor", "Kimura"
                                        # "alignscore" or None
correlation_method = None               # "phylogen"
clustalw_gap_opening = 10               # [0:100]
clustalw_gap_extension = 0.2            # [0:10]
clustalw_distance_matrix = "GONNET"     # "GONNET", "BLOSUM" or "PAM"
muscle_max_iteration = 16               # [2:16]
mafft_configuration = "fftnsi"          # "fftnsi" or "linsi"
alignment_score = "sumofpairs"          # "sumofpairs" or False (* "circularsum")
best_results = 20                       # [1:max(scores)]
results_histogram = True                # True or False
results_heatmap = True                  # True or False
results_structure = "pymol"             # "pymol" or False

