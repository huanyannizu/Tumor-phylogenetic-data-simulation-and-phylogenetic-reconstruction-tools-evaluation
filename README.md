#### To simulate data, run command:
>python Data_Simulation.py repSaveData n_clone n_mutation cell_size_min cell_size_max sample_size read_coverage n_trees

#### Parameters in the above command:
- repSaveData: the folder to save the simulated data
- n_clone: the number of clones in a simulated tree
- n_mutation: the number of mutations in a simulated tree
- cell_size_min: minimum cell population size in a simulated tree
- cell_size_max: maxmum cell population size in a simulated tree
- sample_size: a list of samples sizes 
- read_coverage: a list of read coverages 
- n_trees: the number of simulated trees

For example, to simulate data used in MIPUP (Minimum perfect unmixed phylogenies for multi-sampled tumors via branchings and ILP), run the following command:
>python Data_Simulation.py .../path_to_the_folder_to_save_the_simulated_data  10 100 100 200 [5,10,15,20] [100,1000,10000] 100

#### Output:
the simulated data of each tool will be saved in repositories: Data_MIPUP, Data_LICHeE, Data_AncesTree, Data_CITUP, Data_Treeomics, respectively. In each folder, data is catogorised based on its sample size and read coverage and each catogory is save in a folder named in the format: samplesize_readcoverage. 

In addition, information of each simulated tree is saved in repositoried named in the format: samplesize_readcoverage. In each folder, there are a text file saving clonal information of the tree and a gpickle file saving the tree topology, these information will be used when performing evaluation. 

######################################################################################
### Example of evaluating MIPUP using the 1st dataset of 5 samples and read coverage 100:
<img width="1117" alt="screen shot 2017-10-14 at 7 13 09 pm" src="https://user-images.githubusercontent.com/18735754/31577455-7dbc1f0e-b117-11e7-9e20-4fb761400762.png">
here, 5_100_1_G_Dir.gpickle contains tree topology of the 1st dataset, 5_100_1.txt contains clonal information of the 1st dataset, Data_MIPUP_1_ip_columns.csv and Data_MIPUP_1_ip_tree.dot are the outputs of MIPUP of the 1st dataset.

### Example of evaluating LICHeE using the 1st dataset of 5 samples and read coverage 100:
<img width="1113" alt="screen shot 2017-10-14 at 7 44 54 pm" src="https://user-images.githubusercontent.com/18735754/31577492-38a80116-b118-11e7-8e07-72af3318f2f6.png">
here, Data_LICHeE_1.txt.trees.txt is the output of LICHeE of the 1st dataset.

### Example of evaluating AncesTree using the 1st dataset of 5 samples and read coverage 100:
<img width="1111" alt="screen shot 2017-10-14 at 7 46 38 pm" src="https://user-images.githubusercontent.com/18735754/31577511-7b3f718a-b118-11e7-9490-02e15381a4a7.png">
here, Data_AncesTree.sol and Data_AncesTree.dot are the outputs of AncesTree of the 1st dataset.

### Example of evaluating CITUP using the 1st dataset of 5 samples and read coverage 100:
<img width="1146" alt="screen shot 2017-10-14 at 7 49 04 pm" src="https://user-images.githubusercontent.com/18735754/31577521-c6137954-b118-11e7-83e1-a6e8495b9cd9.png">
here, Results_4_0.txt and node4_tree0.0.dot are the outputs of CITUP of the 1st dataset. Please note CITUP usually has more than one best ouputs, and here Results_4_0 is one of the best results of the 1st dataset.

### Example of evaluating Treeomics using the 1st dataset of 5 samples and read coverage 100:
<img width="1148" alt="screen shot 2017-10-14 at 7 51 27 pm" src="https://user-images.githubusercontent.com/18735754/31577539-2084e68e-b119-11e7-9540-40477f14ccbe.png">
here, Data_5_e=0_01_c0=0_5_af=0_05_mlhtree_1 is the output of Treeomics of the 1st dataset. Please note Treeomics usually has more than one outputs, and here Data_5_e=0_01_c0=0_5_af=0_05_mlhtree_1 is one of the outputs of the 1st dataset.







