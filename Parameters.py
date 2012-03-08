###############################################################################
# Encoding utf-8                                                              #
# Created by F. Madeira, 2012                                                 #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

"Editable parameters"

surface_threshold = 7                   # 7% (Angstrom^2), [0, max(surface)] 
psiblast_evalue = 10                    # [0.0000001:10]
psiblast_identity = 0                   # [0:100] (%)
psiblast_coverage = 0                   # [0:100] (%)
psiblast_threading = False              # Number of cores/servers or False
pairwise_distance = "clustalw"          # "clustalw", "pdistance", "Kimura"
                                        # "jukescantor" or "alignscore"
alignscore_matrix = "BLOSUM62"          # "BLOSUM62" or "PAM250"
theilsen_cutoff = 0.5                   # [0.25:1.0]
clustalw_gap_opening = 10               # [0:100]
clustalw_gap_extension = 0.2            # [0:10]
clustalw_distance_matrix = "GONNET"     # "GONNET", "BLOSUM" or "PAM"
muscle_max_iteration = 16               # [2:16]
mafft_configuration = "linsi"           # "fftnsi" or "linsi"
mafft_threading = False                 # Number of cores/servers or False
alignment_score = "sumofpairs"          # "sumofpairs" or False
best_results = 20                       # [1:max(scores)]
results_histogram = True                # True or False
results_heatmap = True                  # True or False
results_structure = "pymol"             # "pymol" or False

