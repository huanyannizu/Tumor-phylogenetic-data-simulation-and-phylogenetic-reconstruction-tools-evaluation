#### To simulate data, run command:
>python Data_Simulation.py repSaveData n_clone n_mutation cell_size_min cell_size_max n_sample read_cov d

#### Parameters in the above command:
- repSaveData: the directory to save the simulated data
- n_clone: the number of clones(nodes) in a simulated tree
- n_mutation: the number of mutations in a simulated tree
- cell_size_min: minimum cell population size in a simulated tree
- cell_size_max: maxmum cell population size in a simulated tree
- n_sample: the number of samples taken from a simulated tree
- read_cov: read coverage used when sampling from a simulated tree
- d: the number of nodes to remove mutations from

######################################################################################
Example of evaluating MIPUP for dataset "Data_MIPUP_0.txt" (10 nodes,100 mutations,d=0,n_sample = 5,read_coverage = 1000):
![image](https://user-images.githubusercontent.com/18735754/40673985-209b42d8-637c-11e8-938e-602f8a46acd7.png)
Here, 5_1000_0_G_Dir.gpickle contains tree topology of "Data_MIPUP_0.txt", 5_1000_0.txt contains clonal information "Data_MIPUP_0.txt", Data_MIPUP_0_ip_s0_columns.csv and Data_MIPUP_0_ip_s0_tree.dot are the first output of "Data_MIPUP_0.txt".

Example of evaluating LICHeE for dataset "Data_LICHeE_0.txt" (10 nodes,100 mutations,d=0,n_sample = 5, read_coverage = 1000):
![image](https://user-images.githubusercontent.com/18735754/40674534-cf0e75d2-637d-11e8-98c0-47ecf5312fbb.png)
Here, 5_1000_0_G_Dir.gpickle contains tree topology of "Data_LICHeE_0.txt", 5_1000_0.txt contains clonal information "Data_LICHeE_0.txt", Data_LICHeE_0.txt.trees.txt is the output of "Data_LICHeE_0.txt".

Example of evaluating AncesTree for dataset "Data_AncesTree_0.txt" (10 nodes,100 mutations,d=0,n_sample = 5, read_coverage = 1000):
![image](https://user-images.githubusercontent.com/18735754/40674657-22eefd0c-637e-11e8-9653-196e702eafa1.png)
Here, 5_1000_0_G_Dir.gpickle contains tree topology of "Data_AncesTree_0.txt", 5_1000_0.txt contains clonal information "Data_AncesTree_0.txt", Data_AncesTree_0.sol and Data_AncesTree_0.dot are the output of "Data_AncesTree_0.txt".

Example of evaluating CITUP for dataset "Data_CITUP_0.txt" (10 nodes,100 mutations,d=0,n_sample = 5, read_coverage = 1000):
![image](https://user-images.githubusercontent.com/18735754/40674723-511bdefc-637e-11e8-94e5-bf2135b97b2b.png)
Here, 5_1000_0_G_Dir.gpickle contains tree topology of "Data_CITUP_0.txt", 5_1000_0.txt contains clonal information "Data_CITUP_0.txt", Results_5_7.txt and node5_tree7.0.dot are one of the best outputs of "Data_CITUP_0.txt".

Example of evaluating Treeomics for datasets "Data_Treeomics_n_0.txt" and "Data_Treeomics_x_0.txt" (10 nodes,100 mutations,d=0,n_sample = 5, read_coverage = 1000):
![image](https://user-images.githubusercontent.com/18735754/40674879-b6978f7e-637e-11e8-9a56-ad033e117060.png)
Here, 5_1000_0_G_Dir.gpickle contains tree topology of "Data_Treeomics_n_0.txt" and "Data_Treeomics_x_0.txt", 5_1000_0.txt contains clonal information of "Data_Treeomics_n_0.txt" and "Data_Treeomics_x_0.txt", Data_5_e=0_01_c0=0_5_af=0_05_mlhtree_1 is one of the outputs of "Data_Treeomics_n_0.txt" and "Data_Treeomics_x_0.txt".
