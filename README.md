# GWAS

different R function to perform GWAS

The script has been optimized to run with data from Arabidopsis thaliana

gwas.r is a script to run GWAS with a mixed model

The Model uses the fast approximation proposed by Kang et al. 2010
The Null model is solved using the AI algorithm (Gilmour et al. 1995) as it is implemeted in the R library gaston.
 
plots_gwas.r is a script for plotting the results

emma.r use the emma function (Kang et all. 2008) for the MRLE
This is used to fit the full model in the gwas.r script to update the pvalues of the top SNPs.

sim_Y.r is a function to simulate Phenotypes to benchmark and test GWAS approaches. Data need for the simulations are provided in the data section.

ALl data to run GWAS in Arabidopsis can be found at go.uniwue.de/snpextractor
see also aragwas.1001genomes.org

Additional sample data for testing purpose  can be find in the MTMM repository. 
