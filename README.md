# Geometric-renormalization-of-weighted-network
# Codes and data sets for "Geometric renormalization of weighted networks".
Citation: Muhua Zheng, Guillermo García-Pérez, Marián Boguñá, and M. Ángeles Serrano, Geometric renormalization of weighted networks, arXiv:2307.00879 (2023).


In the folder "data", we give the empirical data. There are 3 files for each network:

"**_edgelist.txt"-------edge list file for each line indicates an edge node_i---node_j---weights_wij

"**_coordinates.txt"----coordinates for each node. Each line indicates node kappa theta. Noted that node index has been reasigned by ascending theta. The corresponding edge lists are consistent.  

"**_parameters.txt"-----The parameters:  # N	beta	mu		a		eta		alpha		epsilon2		chi2            

There are 4 C++ codes with different versions of geometric renormalization real network:

(I) 	Renormalization_weight_network_Max.cpp----the code for sup-GRW 

(II)	Renormalization_weight_network_Sum.cpp----the code for sum-GRW 

(III)	Renormalization_weight_network_Phi.cpp----the code for phi-GRW 

(IV)	Renormalization_Binarized.cpp----the code for binarized GR, i.e., it is for the unweighted network. 

(1) Compile the code with "g++ **.cpp -o GR"

(2) Run the code with: ./GR <edgelist> <parameters> <coordinates> <rootname_for_outfile> <total_layers>. edgelist, parameters, and coordinates are the input files for each network from empirical data. <rootname_for_outfile> is the rootname of files for the output. Total_layers is the total layers you want to renormalize. 
 
For example: run the code: ./GR data/cargoshipsBB_edgelist.txt data/cargoshipsBB_Final_parameters.txt data/cargoshipsBB_coordinates.txt data/cargoshipsBB 3
It will renormalize the cargoship network to layer 3.

# For the weighted version(I-III), the codes return some related files with rootname_for_outfile:

rootname_for_outfile+"_C_K_N.txt"-------------------------returns the main statistical properties of each layer. Each line indicates # layer  meanC  meanK  realN  N  ave_weight  ave.strength  mu  
	
rootname_for_outfile+"_RG_layer_**_coordinates.txt"-------is the coordinates and statistics for each node in layer **: # i    kappa    theta    sigma    strength    degree    disparity 
		
rootname_for_outfile+"_RG_layer_**_weight_edgelist.txt"---is the weighed edgelist in layer **. Each line indicates i j wij
		
rootname_for_outfile+"_RG_layer_**_parameters.txt"--------is the parameters in layer **: # N    beta    mu    nu    alpha    eta    a    
		
rootname_for_outfile+"_RG_layer_**_supernode.txt"---------records node i in which supernode.
		
rootname_for_outfile+"_RG_layer_**_Yl.txt"----------------is the disparity of links.
		
# For the binarized version(IV), the code returns some related files with rootname_for_outfile:	

rootname_for_outfile+"_C_K_N.txt"-------------------------returns the main statistical properties of each layer. Each line indicates # layer  meanC  meanK  realN  N  mu  

rootname_for_outfile+"_RG_layer_**_coordinates.txt"-------is the coordinates for each node in layer **: # i    kappa    theta
  
rootname_for_outfile+"_RG_layer_**_edgelist.txt"----------is the edgelist in layer **. Each line indicates i j
  
rootname_for_outfile+"_RG_layer_**_parameters.txt"--------is the parameters in layer **: # N    beta    mu   
  
rootname_for_outfile+"_RG_layer_**_supernode.txt"---------record node i in which supernode.


