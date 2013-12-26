IndEvi and IndEviRe Program. 

This program is used for paper "Transactional Database Transformation and Its Application in Prioritizing Human Disease Genes", Yang Xiang, Philip R.O. Payne, and Kun Huang, IEEE/ACM Transactions on Computational Biology and Bioinformatics 9(1), pp 294-304, 2012. URL: http://dl.acm.org/citation.cfm?id=2077966

Software License Agreement: You may use or modify this computer program for research purposes, provided that you properly cite our paper in publication. This computer program is provided on an as is basis and there is no guarantee on the program nor additional support offered. Neither the author(s) nor their institute(s) is liable under any circumstances. This program archive (including this license agreement) may be updated without further notice.

The following of this README file describes how to run the program.

IndEvi uses the results from "mafia" (licensed under BSD, http://himalaya-tools.sourceforge.net/Mafia/). A binary executable file of MAFIA for Linux is included the in the "bin" directoy for your convenience.

Original Dataset: genes_to_phenotype.txt 
Please get the latest version file from Human Phenotype Ontology Consortium website: http://human-phenotype-ontology.org

!!!Important Notes:
(1) The genes_to_phenotype.txt file, including its content and format, has been updated from time to time by Human Phenotype Ontology Consortium. 
(2) The current release version has passed the test on the genes_to_phenotype.txt downloaded on 12/23/2011.  
	The file format may change such that in the future it is possible the provided Matlab-Preprocessing code cannot handle the file. We are not responsible for maintaining and updating our program.


In the following we illustrate how to use this package.

1. Find and download genes_to_phenotype.txt from Human Phenotype Ontology Consortium website: http://human-phenotype-ontology.org
2. Put genes_to_phenotype.txt in the Matlab-Preprocessing directory.
3. In Matlab-Preprocessing directory using Matlab run preprocessingg2f.m, which will create three preprocessed data files in datasets/ directory 
4. Use Mafia to generate frequent closed itemsets (support level 0.05 is enough to get all closed itemsets):
	bin/mafia -fci 0.0005 -ascii datasets/genes_to_phenotype.dat datasets/freq_itemsets/genes_to_phenotype_0.05_fci.mafia
5. Type "make" to get executable file "IndEvi", and type "make IndEviRe" to get executable file "IndEviRe".
6. Execute bin/IndEvi.  Usage is shown by running it with no parameters.
7. Execute bin/IndEviRe. Usage is shown by running it with no parameters.
	

Example usage for bin/IndEvi (The running time depends on the input file and your machine. For the genes_to_phenotype.txt file release on 12/23/2011, the running time is no more than an hour in our department cluster.):

bin/IndEvi fci datasets/genes_to_phenotype.dat 0.05 787

The last parameter is the query phenotype. 

After executing the command, IndEvi will generate five files:

(1) genes_to_phenotype_support-0.05_fci.max_Cliques
(2) genes_to_phenotype_support-0.05_fci.TopGenes
(3) genes_to_phenotype_support-0.05_fci.score_matrix_covercount
(4) genes_to_phenotype_support-0.05_fci.score_matrix_itmlength
(5) genes_to_phenotype_support-0.05_fci.score_matrix_trwidth

The first file is used by IndEviRe to recover evidence. 
The second file is the gene rank with respect to the query phenotype.
The file (3)-(5) are outputs that can be used (together with genes_to_phenotype.dat) to calculate top genes for any phenotype.

Examples usage for bin/IndEviRe:

bin/IndEviRe fci datasets/genes_to_phenotype.dat 0.05 396 787
bin/IndEviRe fci datasets/genes_to_phenotype.dat 0.05 39 787

If the entry is 0, IndEviRe will output both the maximal biclique (of M) that leads to the maximum F score, and the maximal supporting biclique. The maximal supporting pattern can be easily inferred from the maximal supporting biclique.

If the entry is 1, IndEviRe will output only the maximal biclique (of M) that leads to the maximum F score. This maximal biclique is exactly the maximal supporting pattern (implied by Lemma 3). 
However, to get the maximal supporting biclique you need to remove the query transaction and item (a trivial operation). To avoid essentially duplicate information we do not output it in this version of IndEviRe.


