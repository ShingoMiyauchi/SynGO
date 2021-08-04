# Pair-wise syntenic comparison with SynGO  

                                              Shingo Miyauchi 3Aug21
# Descriptions

Welcome to our visual-omics platform! Synteny Governance Overview (SynGO) is a set of custom R scripts that visualise synteny combined with genomic features from TINGO and PRINGO pipelines (Miyauchi et al., 2020). 
*Please note that SynGO itself does not produce bioinformatic results. It is a visualisation tool for complex genomic data. 

**NOTE 1: Synteny detection**
Identifying syntenic blocks among the species was performed with R package DECIPHER (http://www2.decipher.codes/; see the parameters in Looney et al., 2021). This process is excluded because it is computationally intensive. A high performance computing cluster may be required for calculations. Setting up R environments on the computing cluster is out of scope in this demo.  

**NOTE 2: Genomic features**
We performed prediction of secretome and identification of TEs using in-house pipelines at INRAE Nancy, France (Pellegrin et al. 2015; Morin et al., 2019). Then, output files were cleaned and combined using PRINGO and TINGO pipelines (Miyauchi et al., 2020). 

# Requirement for R and packages

R Studio,
R 3.6.3 

dplyr 1.0.2,
Biostrings 2.54.0,
circlize 0.4.10,
RColorBrewer 1.1-2,
rtracklayer 1.46.0,
wesanderson 0.3.6

# How to run this demo
1) Download INPUT folder 1 to 4 and put all files into one folder 
2) Start up R Studio by clicking an R script
3) Make sure your working directory is in INPUT folder so that R recognises input files. 
4) Read instructions in the scripts and execute the code (line by line).  
5) Figures are generated as output in the same folder. 

# Input files - SynGO_INPUT(1 to 4).zip

1) JGI fungal genomes including scaffold 1 to 10 

JGI_unmasked_assembly_scaffold1to10_fasta
 
2) Coordinates for genes and TEs in scaffold 1 to 10 from TINGO + PRINGO

SynGO_Scaff1to10_GenomeFeature.csv

3) Identified pair-wise synteny among the species with R package DECIPHER

SynGO_synteny_locations
 
# Description of output - Hanabi (fireworks) plots

Hanabi plots show macrosynteny between two fungi with the genomic locations of genes, predicted secretome (e.g. secreted CAZymes, SSPs, lipases, proteases), and TEs in the largest 10 scaffolds. Some selected figures were used in Looney et al., (2021).

**NOTE 3:** Raw figures generated may not be aesthetic enough and it may require adjustments and beautification with Adobe Illustrator.

# References
1. Looney B, Miyauchi S, Morin E, Drula E, Courty PE, Kohler A, Lindquist E, Kuo A, LaButti K, Pangilinan J, et al. 2021. Evolutionary priming and transition to the ectomycorrhizal habit in an iconic lineage of mushroom-forming fungi: is preadaptation a requirement? bioRxiv: 2021.02.23.432530.
2. Miyauchi S, Kiss E, Kuo A, Drula E, Kohler A, Sánchez-García M, Morin E, Andreopoulos B, Barry KW, Bonito G, et al. 2020. Large-scale genome sequencing of mycorrhizal fungi provides insights into the early evolution of symbiotic traits. Nature Communications 11: 1–17.
3. Morin E, Miyauchi S, San Clemente H, Chen EC, Pelin A, de la Providencia I, Ndikumana S, Beaudet D, Hainaut M, Drula E, et al. 2019. Comparative genomics of Rhizophagus irregularis, R. cerebriforme, R. diaphanus and Gigaspora rosea highlights specific genetic features in Glomeromycotina. New Phytologist 222: 1584–1598.
4. Pellegrin C, Morin E, Martin FM, Veneault-Fourrey C. 2015. Comparative analysis of secretomes from ectomycorrhizal fungi with an emphasis on small-secreted proteins. Frontiers in Microbiology 6.
