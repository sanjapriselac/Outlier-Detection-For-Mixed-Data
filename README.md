# Outlier Detection for Mixed-attribute Data

Author: Sanja Priselac 

The code for implementing a selection of available unsupervised outlier detection 
techniques for data with mixed attributes and their comparison in terms of 
effectiveness and efficiency. The analysis is limited to unsupervised scoring techniques where the true status of the
observations is unknown and the output of the method provides scores rather than just
a binary label. The review of scientific literature resulted in eight methods selected
for analysis, designated by the acronyms POD, ABOD, FAMDAD, SECODA, ZDisc,
KMeans, PCAmix, and MIX. Their properties are acquired based on extensive simulation
experiments and evaluation with real data sets. The performance of the methods for
different data structures is investigated by observing the effects of the outlier proportion,
severeness and type, the correlation between attributes, and the different data sizes. 

The implementation of the corresponding methods is based on the following papers:
 * *POD*: "An Effective Pattern Based Outlier Detection Approach for Mixed Attribute Data" (DOI 10.1007/978-3-642-17432-2_13)
 * *ABOD*: "Association-based outlier detection for mixed data” (DOI 10.17485/ijst/2015/v8i25/80260)
 * *FAMDAD*: "Factor analysis of mixed data for anomaly detection" (DOI 10.48550/arXiv.2005.12129)
 * *PCA_T2*: "Outlier Detection using PCA mix based T2 Control Chart for continuous and categorical data" (DOI 10.1080/03610918.2019.1586921)
 * *ZDisc_KMeans*: "Implementation of Numerical Attribute Discretization for Outlier Detection on Mixed Attribute Data Set" (DOI 10.1109/ICOIACT.2018.8350795)
 
The implementation ocde for methods *SECODA* and *MIX* is avaliable on GitHub:
 * *SECODA*: https://github.com/ralfoan/SECODA
 * *MIX*: https://github.com/xuhongzuo/MIX

The code is implemented in the R programming language. 

### Folder structure ###
  * Methods: The implemntation of the outlier detection methods 
  * Simulation.R : The script that generates artificial data for the comparison experiments.
                   Parameters required for data generation include
	* the number of observations n,
	* the total number of attributes p,
	* the number of categorical attributes d,
	* the correlation coefficient ρ,
	* the outlier proportion ϵ,
	* the shift parameter μ,
	* the outlier type t
  * Real_Data.R: The script that performs the comparison of the methods on the selection of real data. All data sets were downloaded from the UCI Machine
	Learning Repository (https://archive.ics.uci.edu/ml/index.php).
	* “Statlog (german credit data) data set.” https://archive.ics.uci.edu/ml/datasets/statlog+(german+credit+data) 
	* “Adult data set.” https://archive.ics.uci.edu/ml/datasets/adult
	* “Heart desease data set.” https://archive.ics.uci.edu/ml/datasets/heart+disease
	* “Contraceptive method choice data set.” https://archive.ics.uci.edu/ml/datasets/Contraceptive+Method+Choice
	* “Thoracic surgery data data set.” https://archive.ics.uci.edu/ml/datasets/Thoracic+Surgery+Data