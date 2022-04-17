# Amino acid network analysis

Comparison of the amino acid networks with different molecular characteristics obtained using molecular dynamics allows us to study the contribution of certain amino acids to both structural stability and conformational mobility of proteins. 

_/Scripts/1_Weight_formula_MD.R_
Considering the entire complex of amino acid parameters in the construction of networks, an algorithm for determining the formula of the node weight has been developed. Weighted graphs were constructed, the weight in which was determined by one of the following parameters –percentage of time amino acid residue was located in forbidden regions, percentage of time amino acid residues interacted with water, RMSF (root mean square fluctuation), conservativeness of the amino acid residue. An unweighted graph was also constructed as a control graph. The parameters of the nodes of the weighted graphs were compared with similar parameters of the nodes of the control unweighted graph using the Kruskal-Wallis test with the pairwise Wilcoxon test adjusted for the multiple testing. The coefficient for each amino acid parameter was calculated based on the significance level obtained during the comparison. The coefficients and values of the parameters determine the weight formula.

_/Scripts/2_Comparative_analysis.R_
To compare amino acid networks adjacency matrices were obtained considering the obtained weights. Clustering using proximity propagation is performed based on adjacency matrices. The similarity score was calculated for each cluster. The similarity score of the graphs was calculated as the median of the similarity scores of clusters. 

_/Scripts/3_DXylose_transporter_comparison.R_
**An example** - The developed algorithm was utilized to conduct a comparative analysis of the amino acid networks of D-xylose transporter of E.coli in different lipid environments.