### Set working directory ###
setwd("~/boxdrive/_Mrug/Diet") 

### Load libraries ###
library(phyloseq)
library(tidyverse)

############################
####### Import data ########
############################

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("./R_analysis/_bin/Diet_asv_table.txt",sep="\t",header=TRUE, row.names=1)
samples<-read.table("./Microbiome_QC/Diet_map_seqreads.txt",sep="\t",header=T,row.names=1)
taxon<-read.table("./R_analysis/_bin/Diet_Silva_tax_table.txt",sep="\t",header=T,row.names=1)
tree <- ape::read.tree("./R_analysis/_bin/Diet.tree")

# prepare phyloseq objects
OTU <- otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX <-tax_table(taxon)
TREE<-phy_tree(tree)
sampledata = sample_data(samples)

# combine
ps <- phyloseq(OTU,TAX,TREE,sampledata)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1253 taxa and 137 samples ]
# sample_data() Sample Data:       [ 137 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 1253 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1253 tips and 1251 internal nodes ]

# mean, max and min of sample read counts
min(sample_sums(ps)) #24341
mean(sample_sums(ps)) #48246.93
max(sample_sums(ps)) #91073
sd(sample_sums(ps))/sqrt(136) #907.2953

# calculate F/B ratio
set.seed(5455)
rare_ps <- phyloseq::rarefy_even_depth(ps)

plot_bar(ps, "Phylum", fill = "Phylum")


phyla <- phyloseq::tax_glom(rare_ps, taxrank = "Phylum")
phyla_rel <- phyloseq::transform_sample_counts(phyla, function(x) { x/sum(x) } )
bac <- as.data.frame(phyloseq::otu_table(phyloseq::subset_taxa(phyla_rel,Phylum=="Bacteroidota")))
firm <-  as.data.frame(phyloseq::otu_table(phyloseq::subset_taxa(phyla_rel,Phylum=="Firmicutes")))
prot <-  as.data.frame(phyloseq::otu_table(phyloseq::subset_taxa(phyla_rel,Phylum=="Proteobacteria")))
verr <-  as.data.frame(phyloseq::otu_table(phyloseq::subset_taxa(phyla_rel,Phylum=="Verrucomicrobiota")))
actin <-  as.data.frame(phyloseq::otu_table(phyloseq::subset_taxa(phyla_rel,Phylum=="Actinobacteriota")))
# F/B Ratio
fb_ratio_log2 <- log2(firm / bac) %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("fb_ratio_log2"= 2)
fb_ratio <- (firm / bac) %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("fb_ratio"= 2)
# Add to sample metadata
sampledata<- full_join(samples,fb_ratio_log2)
sampledata<- full_join(sampledata,fb_ratio)
########
bac <- bac %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("Bacteroidota_ra"= 2)
firm <- firm %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("Firmicutes_ra"= 2)
prot <- prot %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("Proteobacteria_ra"= 2)
verr <- verr %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("Verrucomicrobiota_ra"= 2)
actin <- actin %>% tibble::rownames_to_column(., "Sample_ID") %>% dplyr::rename("Actinobacteriota_ra"= 2)
########
sampledata<- full_join(sampledata,bac)
sampledata<- full_join(sampledata,firm)
sampledata<- full_join(sampledata,prot)
sampledata<- full_join(sampledata,verr)
sampledata<- full_join(sampledata,actin)

row.names(sampledata) <- sampledata$Sample_ID

# export
write.table(sampledata,"~/boxdrive/_Mrug/Diet/R_analysis/FB_ratio/Diet_FB_ratio",sep="\t",col.names=NA, row.names=TRUE)

