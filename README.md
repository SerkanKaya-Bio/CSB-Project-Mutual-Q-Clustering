# CSB-Project-Mutual-Q-Clustering
This repository comprises three scripts to extract C-alpha backbones from protein structures, compute the similarity matrix using mutual-Q, and perform hierarchical agglomerative clustering with centroid linkage on the mutual-Q matrix.

Some RibbonFold-predicted alpha-synuclein 1-99 structures have been added to Ribbons/TestRibbons to run a test analysis.

Here is a step-by-step guide to run the scripts:
1. Put your Ribbon-Structures in the Ribbons folder.
2. In extractCA.py, mutualQ.py and hierarchicalClusteringNew.py, set parameters "protein_name" and "folder" to the correct ones, matching your data.
3. Run scripts in this order and make sure that each script finishes computation before starting the next one (depending on the dataset, running time can be quite long): extractCA.py, mutualQ.py and hierarchicalClusteringNew.py
4. Outputs of extractCA.py, mutualQ.py and hierarchicalClusteringNew.py can be viewed in the folders Calpha, Qmatrix and Clusters, respectively.

Citation:

Guo, L., Yu, Q., Wang, D., Wu, X., Wolynes, P. G., & Chen, M. (2025). Generating the polymorph landscapes of amyloid fibrils using AI: RibbonFold. Proceedings of the National Academy of Sciences, 122(16), e2501321122. https://doi.org/10.1073/pnas.2501321122

Chen, M., Schafer, N. P., & Wolynes, P. G. (2018). Surveying the Energy Landscapes of Aβ Fibril Polymorphism. The Journal of Physical Chemistry. B, 122(49), 11414. https://doi.org/10.1021/acs.jpcb.8b07364

Chen, M., Lin, X., Lu, W., Onuchic, J. N., & Wolynes, P. G. (2016). Protein Folding and Structure Prediction from the Ground Up II: AAWSEM for α/β Proteins. The Journal of Physical Chemistry. B, 121(15), 3473. https://doi.org/10.1021/acs.jpcb.6b09347
