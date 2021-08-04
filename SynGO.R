#  Pair-wise syntenic comparison with SynGO  

#                                   v1.0 Shingo Miyauchi 29July21
#------------------------
#      Descriptions
#------------------------
# <Info>
# Synteny Governance Overview (SynGO) is a set of custom R scripts that visualise synteny combined with genomic features from TINGO and PRINGO pipelines (Miyauchi et al., 2020). Please note that SynGO itself does not produce bioinformatic results. It is a visualisation tool for large genomic data. 

# NOTE 1: Synteny detection
# Identifying syntenic blocks among the species was performed with R package DECIPHER (http://www2.decipher.codes/; see the parameters in Looney et al., 2021). This process is excluded because it is computationally intensive. You may need use a high performance computing cluster for calculations. Regrettably, setting up R environments on the computing cluster is out of scope in this demo.  
# NOTE 2: Genomic features
# We performed prediction of secretome and identification of TEs using in-house pipelines at INRAE Nancy, France (Pellegrin et al. 2015; Morin et al., 2019). Then, clean combined output files were made with PRINGO and TINGO pipelines (Miyauchi et al., 2020). 

# <Input files>
# 1) JGI fungal genomes including scaffold 1 to 10 (unzip file before use!) 
# ###_AssemblyScaffolds.scaffold1to10.fasta

# 2) Coordinates for genes and TEs in scaffold 1 to 10 from TINGO + PRINGO
# SynGO_Scaff1to10_GenomeFeature.csv

# 3) Identified pair-wise synteny among the species with R package DECIPHER
# SynGO_synteny_locations

# <Methods>
# PART 1: Make BED format files for visualisation
# PART 2: Generate Hanabi plots  

# <Requirement for R and packages>
# Please install: 
# R Studio,
# R 3.6.3 
# 
# dplyr 1.0.2
# Biostrings 2.54.0
# circlize 0.4.10
# RColorBrewer 1.1-2
# rtracklayer 1.46.0
# wesanderson 0.3.6
# 
# <How to run the demo>
# 1) Download files into a folder
# 2) Automatically start up R Studio by clicking an R script
# 3) Make sure your working directory is in the folder downloaded so that R recognises the input files. 
# 4) Read instructions in the scripts and execute the code (line by line).  
# 5) Figures are generated as output in the same folder. 
# 
# NOTE 3: Please note that raw figures generated require adjustments and beautification with Adobe Illustator.  
# 
# <Hanabi (fireworks) plots>
# This is the core script integrating the genomic locations of genes and predicted secretome (e.g. secreted CAZymes, SSPs, lipases, proteases), TEs, and synteny between two fungi in the largest 10 scaffolds. Some selected figures were used in Looney et al., (2021).
#------------------------------------
#===============
# Preparation
#===============
# Global option to switch off factorisation. 
options(stringsAsFactors=FALSE)

# Load packages
library(dplyr)
library(Biostrings)
library(circlize)
library(RColorBrewer)
library(rtracklayer)
library(wesanderson)

# Colour functions
colour_fun <- function(n, color.set) {
    col.ls<-colorRampPalette(brewer.pal(brewer.pal.info [which(rownames(brewer.pal.info)==color.set),"maxcolors"],color.set))(n)
    return(col.ls[1:n])}

############################################
# PART 1: Make BED format inputs for Circos
############################################
#======================================
# 1) Get pair-wise comparison of fungi
#======================================
# Load synteny information 
syn.df <- read.csv(list.files(pattern="SynGO_synteny_locations"))

# Get pair-wise comparison IDs
syn.df$pair <- paste0(syn.df$Fungus.1, ".", syn.df$Fungus.2)

# Get unique pair-wise IDs 
pair.id <- unique(syn.df$pair)

#=======================
# 2) Get scaffold 1 to 10 
#=======================
# Get scaffold IDs from fasta files
fasta.fl <- list.files(pattern="*AssemblyScaffolds.scaffold1to10.fasta",recursive =T) %>% print() 

# Get fungal ID from fasta files
fungi.ID <- unique(gsub("_AssemblyScaffolds.scaffold1to10.fasta",
                        "", basename(fasta.fl))) %>% print()

# Read and make a data frame with scaffold 1 - 10
scaf.bed.df <- data.frame()
for (j in 1:length(fungi.ID)) {
    # Loading fasta scaffold 1 - 10
    fast <- readDNAStringSet(fasta.fl[j])
    # list of data frames 
    df <- data.frame(ID= fungi.ID[j], 
                     scaffold = names(fast), 
                     start = 1,
                     end = width(fast))
    scaf.bed.df <-rbind(scaf.bed.df, df)
}

# Make a new column with fungal ID + scaffold # 
scaf.bed.df$scf <- paste0(substring(gsub("\\d", "", scaf.bed.df$ID) , 1, 4), 
                          gsub("\\w+_|\\w+_0", "", scaf.bed.df$scaffold))

#--------------------------------------
# Make a pair-wise bed list for scaffolds
#--------------------------------------
scaf.circos.ls <- list()
for (k in 1:length(pair.id)) {
    # Extract pair-wise Fungal ID  
    fungus1 <- unlist(strsplit(pair.id[k], "\\."))[1]
    fungus2 <- unlist(strsplit(pair.id[k], "\\."))[2]
    
    # Make a bed list
    scaf.circos.ls[[k]] <- rbind(scaf.bed.df[scaf.bed.df$ID == fungus1, ], 
                                 scaf.bed.df[scaf.bed.df$ID == fungus2, ])
}

# Transfer pair.id to the list
names(scaf.circos.ls) <- pair.id

# Remove temporary objects
rm(df, j, fasta.fl, fast, fungi.ID, k)

#===============================
# 3) Genes in scaffold 1 to 10 
#===============================
# Load a table with coordinates of genome features  
gene.df <- read.csv(list.files(pattern="SynGO_Scaff1to10_GenomeFeature.csv"))

# Subset by protein type - cazy, ssp, lip.pro, non_secretome
gene.bed.df <- gene.df[gene.df$type %in% c("cazy","SSP", "lip.pro", "non_secretome"),]

#---------------------
# Calculate distances between the genes
#---------------------
# Get a list of single fungus
ID.ls <- unique(scaf.bed.df$ID)

# Calculate distances per fungus
gene.sub.df.ls <- list()
for (j in 1:length(ID.ls)) {
    # Subset gene.bed.df by a fungus
    gene.sub.df <- gene.bed.df[gene.bed.df$fungus == ID.ls[j],]
    
    # Calculate distances between all TEs
    gene.sub.dist <- rainfallTransform(gene.sub.df[, c("scf", "start", "end")], mode="mean")
    
    # IMPORTANT: It fails to get numbers some times because of the no distance caused by having two gene names like SSP + CAZyme - in that case add 1
    if (any(is.na(gene.sub.dist$dist))) {gene.sub.dist$dist[is.na(gene.sub.dist$dist)] <- 1}
    
    # Transform into log10
    gene.sub.dist$dist <- log10(gene.sub.dist$dist +1)
    
    # Column bind distances back onto gene.sub.df
    gene.sub.df <- cbind(gene.sub.df, dist= gene.sub.dist$dist)
    
    # list up all data frames
    gene.sub.df.ls[[j]] <- gene.sub.df
}

# Row bind all data frames
gene.bed.df <- do.call(rbind, gene.sub.df.ls)

# Get max distance for genes for limit of Y axis in Hanabi plot
gene.dist.max <- max(gene.bed.df$dist)

#--------------------
# Separate gene.bed.df into CAZy, SSP, Lipase, Protease, non-secreted
#---------------------
# Extract bed for CAZy, SSPs, lipases+proteases,  non-secretome, and TEs
CAZy.bed.df <- gene.bed.df[gene.bed.df$type == "cazy",]
SSP.bed.df <- gene.bed.df[gene.bed.df$type == "SSP",]
LipPro.bed.df <- gene.bed.df[gene.bed.df$type == "lip.pro",]
Non_secretome.bed.df <- gene.bed.df[gene.bed.df$type == "non_secretome",]

# Remove temporary objects
rm(gene.sub.df.ls, gene.bed.df, gene.sub.dist, j)

#--------------------------------------
# Make a pair-wise bed list for the genes
#--------------------------------------
gene.circos.ls <- list()
for (k in 1:length(pair.id)) {
    # Extract pair-wise Fungal ID  
    fungus1 <- unlist(strsplit(pair.id[k], "\\."))[1]
    fungus2 <- unlist(strsplit(pair.id[k], "\\."))[2]
    
    # Make a bed list
    gene.circos.ls[[k]] <- list(
        cazy = rbind(CAZy.bed.df[CAZy.bed.df$fungus == fungus1,], 
                     CAZy.bed.df[CAZy.bed.df$fungus == fungus2,]),
        
        lipro = rbind(LipPro.bed.df[LipPro.bed.df$fungus == fungus1,], 
                      LipPro.bed.df[LipPro.bed.df$fungus == fungus2,]),
        
        ssp = rbind(SSP.bed.df[SSP.bed.df$fungus == fungus1,], 
                    SSP.bed.df[SSP.bed.df$fungus == fungus2,]),
        
        non_secre= rbind(Non_secretome.bed.df[
            Non_secretome.bed.df$fungus == fungus1,], 
            Non_secretome.bed.df[Non_secretome.bed.df$fungus == fungus2,]))
}

# Name the list of genes 
names(gene.circos.ls) <- pair.id

# Remove intermediate objects
rm(CAZy.bed.df, LipPro.bed.df, Non_secretome.bed.df, SSP.bed.df, k)

#=================================================
# 4) Transposable elements + unknown repeat sequences
#=================================================
# Extract TEs and unkown repeats from gene.df
TE.bed.df <- gene.df[gene.df$type == "TE/Repeat sequences",]

#---------------------
# IMPORTANT: Calculate distances between TEs/unknown repeats
#---------------------
TE.sub.df.ls <- list()
for (j in 1:length(ID.ls)) {
    
    # Subset TE.bed.df by a fungus
    TE.sub.df <- TE.bed.df[TE.bed.df$fungus == ID.ls[j],]
    
    # Calculate distances between all TEs
    TE.sub.dist <- rainfallTransform(TE.sub.df[, c("scf", "start", "end")], mode="mean")
    
    # IMPORTANT: Somehow it fails to get numbers some times - in that case add 1
    if (any(is.na(TE.sub.dist$dist))) 
    {TE.sub.dist$dist[is.na(TE.sub.dist$dist)] <- 1}
    
    # Transform into log10
    TE.sub.dist$dist <- log10(TE.sub.dist$dist +1)
    
    # Row bind distances to TE.sub.df
    TE.sub.df <- cbind(TE.sub.df, dist= TE.sub.dist$dist)
    
    # Make a list of data frames 
    TE.sub.df.ls[[j]] <- TE.sub.df
}

# Put them together to make TE.bed.df
TE.bed.df <- do.call(rbind, TE.sub.df.ls)

# Get max distance for genes for limit of Y axis in Hanabi plots
TE.dist.max <- max(TE.bed.df$dist)

# Remove temporary objects
rm(TE.sub.df.ls, TE.sub.df, TE.sub.dist, j)

#=================================================
# 5) Make a pair-wise bed list for TEs in scaffold 1 to 10
#=================================================
TE.circos.ls <- list()
for (k in 1:length(pair.id)) {
    # Extract pair-wise Fungal ID  
    fungus1 <- unlist(strsplit(pair.id[k], "\\."))[1]
    fungus2 <- unlist(strsplit(pair.id[k], "\\."))[2]
    
    # Make a bed list
    TE.circos.ls[[k]] <- rbind(TE.bed.df[TE.bed.df$fungus == fungus1 & 
                                   TE.bed.df$scf %in% paste0(substring(fungus1, 1, 4), 
                                                           seq(1:10)),], 
                               TE.bed.df[TE.bed.df$fungus == fungus2 & 
                                   TE.bed.df$scf %in% paste0(substring(fungus2, 1, 4), 
                                                           seq(1:10)),])
    }


# Name the list of genes 
names(TE.circos.ls) <- pair.id

# Remove temporary objects
rm(TE.bed.df, k)

#===================================
# 6) Make bed objects for Syntenic regions 
#===================================
# Make a scf column
syn.df$scf.1 <- paste0(substring(syn.df$Fungus.1, 1,4), syn.df$Scaffold.1)
syn.df$scf.2 <- paste0(substring(syn.df$Fungus.2, 1,4), syn.df$Scaffold.2)

# Make a pari-wise list for syntenic regions 
syn.circos.ls <- list()
for (k in 1:length(pair.id)) {
        # Subset by pair wise comparison 
        df <- syn.df[syn.df$pair == pair.id[k],]
        
        # Extract pair-wise Fungal ID  
        fungus1 <- unlist(strsplit(pair.id[k], "\\."))[1]
        fungus2 <- unlist(strsplit(pair.id[k], "\\."))[2]
    
        # Convert pair-wise comparison into bed format for Species 1  
        syn.bed1 = df[df$Fungus.1 == fungus1, 
                      c("scf.1", "Start.1", "End.1")]
        # Same for Species 2
        syn.bed2 = df[df$Fungus.2 == fungus2, 
                      c("scf.2", "Start.2", "End.2")]
        
        # Change the names 
        names(syn.bed1) <- c("scf", "start", "end")
        names(syn.bed2) <- c("scf", "start", "end")

        # Make a bed list
        syn.circos.ls[[k]] <- list(syn.bed1=syn.bed1, syn.bed2 = syn.bed2)
        }

# Name the list 
names(syn.circos.ls) <- pair.id 

# Remove temporary objects
rm(df, syn.bed1, syn.bed2, syn.df, fungus1, fungus2, k)

#################################
# PART 2: Generate Hanabi Plots  
#################################
# NOTE: This process takes a few minutes. It generates 138 combinations of "fungus vs fungus" Hanabi plots. 

for (i in 1:length(pair.id)) {

    # Make bed files per pair-wise comparison    
    #----------------------
    #(1) Scaffold 1 to 10 
    #----------------------
    scaf.bed <- scaf.circos.ls[[i]][c("scf", "start", "end")] 
    
    #-----------------------------------    
    #(2) Secretome + non-secretome in scaffold 1 to 10
    #-----------------------------------
    CAZy.bed <- gene.circos.ls[[i]]$cazy[c("scf", "start", "end", "dist")]
    SSP.bed <-  gene.circos.ls[[i]]$ssp[c("scf", "start", "end", "dist")]
    LipPro.bed <- gene.circos.ls[[i]]$lipro[c("scf", "start", "end", "dist")]
    Non_secretome.bed <- gene.circos.ls[[i]]$non_secre[c("scf", "start", "end", "dist")]
    
    # Make a list CAZy, SSP, LipPro, non-secretome 
    gene.bed.ls = list(Non_secretome.bed, CAZy.bed, SSP.bed, LipPro.bed)
    
    # Name the list 
    names(gene.bed.ls) <- c("NonSecreted", "CAZy", "SSP", "LipPro")

    #------------------------------------
    #(3) TE families + unknown repeat sequences in scaffold 1 to 10 
    #------------------------------------
    #--------------------
    # Separate TE.bed.df into TE families
    #---------------------
    # Get unique TEs per pair-wise comparison
    TE.ls <- unique(TE.circos.ls[[i]]$annot)
    
    # Sort TE.ls - to be consistent for compariosns 
    TE.ls <- sort(TE.ls)
    
    # Index individual TEs per fungus
    TE.bed.ls <- list()
    for (k in 1:length(TE.ls)) {
        TE.bed.ls[[k]] <- TE.circos.ls[[i]][TE.circos.ls[[i]]$annot == TE.ls[k], 
                                    c("scf", "start", "end", "dist")]
    }
    
    # Name the list 
    names(TE.bed.ls) <- TE.ls
    
    #------------------
    # Change the order of Unknown repeat seqs
    #------------------
    # Separate Unknown repeat 
    Unknown.bed <- TE.bed.ls[names(TE.bed.ls) == "Unknown"]
    
    # Remove Unknown from TE.bed.ls and put it back
    TE.bed.ls <- TE.bed.ls[names(TE.bed.ls) != "Unknown"]
    TE.bed.ls <- c(Unknown.bed, TE.bed.ls)
    
    #--------------------
    # bed for pair-wise syntenic regions
    #--------------------
    syn.bed1 <- syn.circos.ls[[i]]$syn.bed1
    syn.bed2 <- syn.circos.ls[[i]]$syn.bed2
    
    #=======================
    # Make + Save Hanabi plots
    #=======================
    # Automatic saving function in PDF
    pdf (file=paste0(getwd(), "/" ,  pair.id[i], "_Synteny_Genes_TEs_", "Scaffold_1to10", ".pdf"), 
         height = 7, width = 7)
    
    # Colours for non-secreted gene (grey95) and secretomes (colours)
    gene.color.ls <- add_transparency(c("grey95", "gold", "limegreen", "blue"), 0.2)
    
    # Get fashonable colours for TEs 
    # NOTE: unknown has a manual colour so -1
    col.num <- length(TE.ls)-1
    TE.color.ls <- add_transparency(c("cornsilk3", wes_palette("Darjeeling1", col.num, type = "continuous")), 0.2)
    
    #-------------------------
    # Make a circle of labells 
    #-------------------------
    circos.genomicInitialize(scaf.bed, plotType = c("axis","labels"))
    
    #-------------------------
    # Make a first track
    #-------------------------
    circos.track(ylim = c(0, 0.5),
                 # Make two identical grey gradients for 10 for fungal species  
                 bg.col = c(colour_fun(20, "Greys")[12:3],
                            colour_fun(20, "Greys")[12:3]), 
                 bg.border = NA, track.height = 0.05)
    
    #-------------------------
    # Second track - CAZy, SSP, lipase, protease, non-secreted genes
    #-------------------------
    circos.genomicTrackPlotRegion(gene.bed.ls,
                                  panel.fun = function(region, value, ...) {
                                      i = getI(...)
                                      circos.genomicPoints(region, value, 
                                                           pch = 20, 
                                                           cex = 0.8, 
                                                           col = gene.color.ls[i],...)},
                                  ylim = c(0, gene.dist.max), 
                                  bg.border = NA, 
                                  track.height = 0.1)
    #-------------------------
    # Third track - TE and unkonwn repeats
    #-------------------------
    circos.genomicTrackPlotRegion(TE.bed.ls,
                                  panel.fun = function(region, value, ...) {
                                      i = getI(...)
                                      circos.genomicPoints(region, value, 
                                                           pch = 20, 
                                                           cex = 0.8, 
                                                           col = TE.color.ls[i],...)},
                                  ylim = c(0, TE.dist.max), 
                                  bg.border = NA, track.height = 0.15)
    #---------------------------
    # Links for syntenic regions 
    #---------------------------
    circos.genomicLink(syn.bed1, syn.bed2, 
                       # Make colours for ribbons
                       col = add_transparency("grey65", 0.4),
                       border = NA)

    # Turn off Circos 
    circos.clear()
    
    #----------------------
    # Add title in the middle
    #----------------------
    # Extract pair-wise Fungal ID  
    fungus1 <- unlist(strsplit(pair.id[i], "\\."))[1]
    fungus2 <- unlist(strsplit(pair.id[i], "\\."))[2]
    
    text(0, -0.3, paste0(fungus1, "vs", fungus2), cex = 1.8, col= add_transparency("grey60", 0.2))
    
    # Add legend 1
    legend(0, 0.5, legend= c(names(gene.bed.ls)),
           col=gene.color.ls, 
           pch= 20, cex=0.75,
           bty = "n")

    # Add legend 2
    legend(-0.4, 0.5 , legend= c(names(TE.bed.ls)[1:(length(names(TE.bed.ls)))]),
           col=TE.color.ls, 
           pch= 20, cex=0.75,
           bty = "n")

    dev.off()
}
