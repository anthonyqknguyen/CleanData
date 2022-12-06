RNA Sequencing Power Analysis


This script can be sourced using the function source("PathTo/RNASeqPower.R") and its functions utilized to perform 
RNA Seq Power Analysis. Requires the installation of the following packages: edgeR, MASS, and dplyr



Functions:


ofPowerAnalysis(counts, factor, sims = 5, nmin = 2, nmax = 10, interval = 1)

	Performs power analysis of a one factor experiment using the provided counts table and factor. 
	Calculates power values for various number of samples ranging from nmin to nmax with the given interval. 
	The power calculation is repeated sims time to obtain an average power for each distinct number of samples.

	Parameters:
		counts - A counts table of RNA Seq data. Columns represent individual samples while rows represent genes.
		factor - The factor with which to perform the experiment with. This factor should take the form as a list
			of discrete number or character values, with each unique character representing a different value
			of the factor. The index of these characters must align with the columns of the counts table.
			Ex. factor = c("A", "A", "A", "B", "B", "B") where A represents smokers and B represents non-smokers
			    in this example, the first 3 columns of the counts table are smokers and the last three columns are
			    non-smokers
		sims - The number of simulations to be performed for each distinct number of samples. This allows for an average
			power value to be obtained for the different numbers of samples.
		nmin - The minimum value of samples with which to calculate the power value. Must be at least 2.
		nmax - The maximum value of samples with which to calculate the power value.
		interval - The interval between numbers of samples with which to calculate power. Increasing the interval will
			improve runtime.
			Ex. with nmin = 5, nmax = 20, and interval = 5, the power will be calculated for 5, 10, 15, and 20 samples.

	Returns:
		results - A matrix with sims columns where each row represents a distinct number of samples as specified in the
			parameters. The entries within the matrix contain the power values.



pairedPowerAnalysis(counts, factor1, factor2, sims = 5, nmin = 2, nmax = 10, interval = 1)

	Performs power analysis of a two factor experiment using the provided counts table and factors. 
	Calculates power values for various number of samples ranging from nmin to nmax with the given interval. 
	The power calculation is repeated sims time to obtain an average power for each distinct number of samples.

	Parameters:
		counts - A counts table of RNA Seq data. Columns represent individual samples while rows represent genes.
		factor1 - The first factor with which to perform the experiment with. This factor should take the form as a list
			of discrete number or character values, with each unique character representing a different value
			of the factor. The index of these characters must align with the columns of the counts table.
			Ex. factor = c("A", "A", "A", "B", "B", "B") where A represents smokers and B represents non-smokers
			    in this example, the first 3 columns of the counts table are smokers and the last three columns are
			    non-smokers
		factor2 - The second factor with which to perform the experiment with. This factor should have the same format as
			factor1
		sims - The number of simulations to be performed for each distinct number of samples. This allows for an average
			power value to be obtained for the different numbers of samples.
		nmin - The minimum value of samples with which to calculate the power value. Must be at least 2.
		nmax - The maximum value of samples with which to calculate the power value.
		interval - The interval between numbers of samples with which to calculate power. Increasing the interval will
			improve runtime.
			Ex. with nmin = 5, nmax = 20, and interval = 5, the power will be calculated for 5, 10, 15, and 20 samples.

	Returns:
		results - A matrix with sims columns where each row represents a distinct number of samples as specified in the
			parameters. The entries within the matrix contain the power values.



estimateOFParams(counts, factor)

	Uses the edgeR package to estimate the negative binomial parameters of the counts table using the given factor. Assumes
	The counts table has a negative binomial distribution.


	Parameters:
		counts - A counts table of RNA Seq data. Columns represent individual samples while rows represent genes.
		factor - The factor with which to estimate the parameters with. This factor should take the form as a list
			of discrete number or character values, with each unique character representing a different value
			of the factor. The index of these characters must align with the columns of the counts table.
			Ex. factor = c("A", "A", "A", "B", "B", "B") where A represents smokers and B represents non-smokers
			    in this example, the first 3 columns of the counts table are smokers and the last three columns are
			    non-smokers

	Returns:
		params - A list where each entry is another list containing parameters for each gene. The entries of params 
			are as follows:
			params[1] - returns y, which is created with the DGEList edgeR function using the counts table and factor.
			params[2] - returns fc, a list of fold change values for each gene
			params[3] - returns dispsCR, a list of dispersion values for each gene
			params[4] - returns sample_data, a data frame containing the factor type as well as the average library
				size for each gene
			params[5] - returns nofit,  a list of indices of genes which the edgeR package could not fit parameters to
			params[6] - returns de, a list of the genes which were found to be differentially expressed in the experiment



estimatePairedParams(counts, factor1, factor2)

	Uses the edgeR package to estimate the negative binomial parameters of the counts table using the given factors. Assumes
	The counts table has a negative binomial distribution.


	Parameters:
		counts - A counts table of RNA Seq data. Columns represent individual samples while rows represent genes.
		factor1 - The first factor with which to perform the experiment with. This factor should take the form as a list
			of discrete number or character values, with each unique character representing a different value
			of the factor. The index of these characters must align with the columns of the counts table.
			Ex. factor = c("A", "A", "A", "B", "B", "B") where A represents smokers and B represents non-smokers
			    in this example, the first 3 columns of the counts table are smokers and the last three columns are
			    non-smokers
		factor2 - The second factor with which to perform the experiment with. This factor should have the same format as
			factor1

	Returns:
		params - A list where each entry is another list containing parameters for each gene. The entries of params 
			are as follows:
			params[1] - returns y, which is created with the DGEList edgeR function using the counts table and factor.
			params[2] - returns fc, a list of fold change values for each gene
			params[3] - returns dispsCR, a list of dispersion values for each gene
			params[4] - returns sample_data, a data frame containing the factor type as well as the average library
				size for each gene
			params[5] - returns nofit,  a list of indices of genes which the edgeR package could not fit parameters to
			params[6] - returns de, a list of the genes which were found to be differentially expressed in the experiment



simulateOF(params, n)

	Simulated negative binomial (counts) data using the given parameters. Simulates 2n samples, n samples representing the first
	value of the factor used to create the parameters, and n samples representing the second value of the factor.

	Parameters:
		params - The object returned by the estimateOFParams function. See the documentation on how to use the function.
		n - The number of samples to be simulated

	Returns:
		m - A data frame of simulated data. The first n columns represent simulated samples of the first value of the factor,
			while the last n columns represent simulated samples of the second value of the factor.



simulatePaired(params, n)

	Simulated negative binomial (counts) data using the given parameters. Simulates 2n samples with columns alternating between
	both factors.

	Parameters:
		params - The object returned by the estimatePairedParams function. See the documentation on how to use the function.
		n - The number of samples to be simulated

	Returns:
		m - A data frame of simulated data. Simulates 2n samples with columns alternating between both factors.
			For example, if the first factor has values "A" and "B", while the second factor has values "1" and "2",
			then the factor of the columns can be represented as "1A", "1B", "2A", "2B" and so on for 2n samples.



evalOFData(m, n)

	Evaluates simulated one factor counts data in order to find which genes are differentially expressed.

	Parameters:
		m - The simulated data returned by the simulateOF function. See the documentation on how to use the function.
		n - The number of samples. This value should be the same as the n value used when simulating the data.

	Returns: 
		pval_list_sim - A list of p-values for each gene, showing the confidence that each gene is differentially expressed.



evalPairedData(m, n)

	Evaluates simulated two factor counts data in order to find which genes are differentially expressed.

	Parameters:
		m - The simulated data returned by the simulatePaired function. See the documentation on how to use the function.
		n - The number of samples. This value should be the same as the n value used when simulating the data.

	Returns: 
		pval_list_sim - A list of p-values for each gene, showing the confidence that each gene is differentially expressed.






