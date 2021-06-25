# ECSWalk
ECSWalk: A carcinogenic driver module detection method based on a network model.

##1. Data used in this study:
In our research, we used data processed by Rafsan Ahmed, et al. and Mark D M Leiserson, et al. Download their research here respectively: https://doi.org/10.1093/bioinformatics/btz655 and https://doi.org/10.1038/ng.3168.

##2. Input：
hint_edge_file.txt: The PPI network file.
hint_index_file.txt: The nodes to index file mapping.
genes_vs_patients_indices_gen_paran.txt: Gene and its mutation patients.
patients_to_indices_gen.txt: The patients name and its patient ID.
pan12gene2freq.txt: The gene mutation frequency.
Census_allTue_May_23_12-08-15_2017.tsv: The cosmic genes.

##3. Run:
compute_edge_weights.py: Compute the weights of edges between gene nodes in biological networks.
random_walk.py: Perform a random walk on a directed weighted network.
connected_components_isolarge.py: Driver module detection.

##4. Outputs:
deg_ed.txt: Mutual exclusion between connected gene nodes.
deg_cd.txt: Coverage between connected gene nodes.
deg_sim.txt: Similarity between connected gene nodes.
ecswalk.txt: Weights of edges between connected gene nodes in biological networks.
ECSWalk folder: Driver module set at different total_size thresholds.
