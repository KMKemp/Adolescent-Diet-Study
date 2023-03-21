
### Set working directory ###
setwd("~/boxdrive/_Mrug/Diet") #desktop
setwd("~/Box/_Mrug/Diet") #laptop

### Load libraries ###
library(tidyverse)
library(phyloseq)
library(mirlyn)
library(btools)


####################################################################
################### Import into Phyloseq for QC ####################
####################################################################
# set files paths
otu <- read.table("./R_analysis/_bin/Diet_asv_table.txt",sep="\t",header=TRUE, row.names=1)
otu <- t(otu)
samples<-read.table("./R_analysis/_bin/Diet_map.txt",sep="\t",header=T,row.names=1)
tree <- ape::read.tree("./R_analysis/_bin/Diet.tree")

# make phyloseq object
Diet <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), phy_tree(tree),
                  sampledata = sample_data(samples))
Diet
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1253 taxa and 137 samples ]
# sample_data() Sample Data:       [ 137 samples by 19 sample variables ]
# phy_tree()    Phylogenetic Tree: [ 1253 tips and 1251 internal nodes ]

####################################################################
######### Calculate alpha diversity and add to mapping file ########
####################################################################

min(sample_sums(Diet2))
# [1] 24341

# Use mirlyn to rarefy
Rarefy_whole_rep_example <- mirlyn::rarefy_whole_rep(Diet2, rep = 100, mc.cores=8)
# View rarefaction curve
Rarecurve_ex <- mirlyn::rarecurve(Rarefy_whole_rep_example, sample = "Sample")
Rarecurve_ex
# set rarefaction to default minimum library size with 1000 iterations and no replacement
diet_mirl <- mirlyn::mirl(Diet2, rep = 100, set.seed = NULL, trimOTUs = FALSE, replace = FALSE, mc.cores = 8)
##### SAVE THE FILE SO YOU DON'T HAVE TO REPEAT
saveRDS(diet_mirl, file="~/boxdrive/_Mrug/Diet/R_analysis/_bin/DietW1_diet_mirl.rds")
diet_mirl<-readRDS("~/boxdrive/_Mrug/Diet/R_analysis/_bin/DietW1_diet_mirl.rds")

shannon <- alphadivDF(diet_mirl, diversity = "shannon")
simpson <- alphadivDF(diet_mirl, diversity = "simpson")
invsimpson <- alphadivDF(diet_mirl, diversity = "invsimpson")

shannon_avg <-Rmisc::summarySE(shannon, measurevar="DiversityIndex", groupvars=c("Sample_ID"))
simpson_avg <-Rmisc::summarySE(simpson, measurevar="DiversityIndex", groupvars=c("Sample_ID"))
invsimpson_avg <-Rmisc::summarySE(invsimpson, measurevar="DiversityIndex", groupvars=c("Sample_ID"))

# Faith's Phylogenetic Diversity (SR = observed species richness)
out <- lapply(diet_mirl, btools::estimate_pd)
PD = (mapply(function(x) cbind(x[1]), out, USE.NAMES = TRUE, SIMPLIFY = FALSE))
SR = (mapply(function(x) cbind(x[2]), out, USE.NAMES = TRUE, SIMPLIFY = FALSE))

PD <-list.cbind(PD) 
SR <-list.cbind(SR)

new_names<-1:1000
colnames(SR) = paste(colnames(SR),new_names)
colnames(PD) = paste(colnames(PD),new_names)

PD$Sample_ID <- rownames(PD)
SR$Sample_ID <- rownames(SR)

SR <- SR %>% mutate(richness = rowMeans(across(where(is.numeric))))
PD <- PD %>% mutate(faithspd = rowMeans(across(where(is.numeric))))

Mirl_alpha <-as.data.frame(cbind(Sample_ID = shannon_avg$Sample_ID,m.shannon = shannon_avg$DiversityIndex, 
                                 m.simpson = simpson_avg$DiversityIndex, m.invsimpson = invsimpson_avg$DiversityIndex, 
                                 faith.pd = PD$faithspd, richness = SR$richness))

alpha <- full_join(Bw_rich, Mirl_alpha, by ="Sample_ID")
write.csv(alpha, "./R_analysis/Alpha_diversity/W1_alphadiv_bw_mirl.csv" )




