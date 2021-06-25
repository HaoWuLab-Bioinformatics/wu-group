# ECSWalk
ECSWalk: A carcinogenic driver module detection method based on a network model.

## **Data used in this study:**

In our research, we used data processed by Rafsan Ahmed, et al. and Mark D M Leiserson, et al. Download their research here respectively: https://doi.org/10.1093/bioinformatics/btz655 and https://doi.org/10.1038/ng.3168.

## **Input：**

1. hint_edge_file.txt: The PPI network file.
2. hint_index_file.txt: The nodes to index file mapping.
3. genes_vs_patients_indices_gen_paran.txt: Gene and its mutation patients.
4. patients_to_indices_gen.txt: The patients name and its patient ID.
5. pan12gene2freq.txt: The gene mutation frequency.
6. Census_allTue_May_23_12-08-15_2017.tsv: The cosmic genes.

## **Run:**
1. compute_edge_weights.py: Compute the weights of edges between gene nodes in biological networks.
2. random_walk.py: Perform a random walk on a directed weighted network.
3. connected_components_isolarge.py: Driver module detection.

## **Outputs:**

1. deg_ed.txt: Mutual exclusion between connected gene nodes.
2. deg_cd.txt: Coverage between connected gene nodes.
3. deg_sim.txt: Similarity between connected gene nodes.
4. ecswalk.txt: Weights of edges between connected gene nodes in biological networks.
5. ECSWalk folder: Driver module set at different total_size thresholds.

