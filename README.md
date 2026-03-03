# CSB-Project-Mutual-Q-Clustering
This repository comprises three scripts to extract C-alpha backbones from protein structures, compute the similarity matrix using mutual-Q, and perform hierarchical agglomerative clustering with centroid linkage on the mutual-Q matrix.

Some RibbonFold-predicted alpha-synuclein 1-99 structures have been added to Ribbons/TestRibbons to run a test analysis.

Here is a step-by-step guide to run the scripts:
1. Put your Ribbon-Structures in the Ribbons folder.
2. In extractCA.py, mutualQ.py and hierarchicalClusteringNew.py, set parameters "protein_name" and "folder" to the correct ones, matching your data.
3. Run scripts in this order and make sure that each script finishes computation before starting the next one (depending on the dataset, running time can be quite long): extractCA.py, mutualQ.py and hierarchicalClusteringNew.py
4. Outputs of extractCA.py, mutualQ.py and hierarchicalClusteringNew.py can be viewed in the folders Calpha, Qmatrix and Clusters, respectively.
