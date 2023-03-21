### Raw Data ###

# The original sequencing reads are available for download from: 
# <https://data.genome.uab.edu/microbiome/result2021/M274-Mrugrerun1/M274_analysis/MrugW1_analysis/ANALYSIS/microbiome_report.html>.

### Set working directory ###
setwd("~/boxdrive/_Mrug/Diet") #desktop
setwd("~/Box/_Mrug/Diet") #laptop

### Load libraries ###

library(dada2)
library(phyloseq)
library(tidyverse)
library(DECIPHER)
library(phangorn)

###############################
########### Purpose ###########
###############################

# The purpose of this script is to re-assign taxo nomy to the W1 diet study microbiome 
# using the latest SILVA and RDP databases because the taxonomy assigned by the core 
# is a #%!$ing mess!!!

#####  Load data ##### 
# this is the raw seqtab.nochim.csv produced by the core and downloaded from the weblink above
seqtab.nochim <- read.csv("~/boxdrive/_Mrug/Diet/Microbiome_QC/raw/seqtab.nochim.csv") 
seqtab.nochim <- data.frame(seqtab.nochim, row.names = 1)
seqtab.nochim <- as.matrix (seqtab.nochim)

#####  Plot the sequence variant count per sequence length
table <- as.data.frame(table(nchar(colnames(seqtab.nochim))))
colnames(table) <- c("LENGTH","COUNT")
ggplot(table,aes(x=LENGTH,y=COUNT)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8)) +
  theme(axis.text.y=element_text(size=10))

##### Remove columns with lengths outside the expected sequence length window
seqlens  <- nchar(colnames(seqtab.nochim))
seqtab.filt <- seqtab.nochim[, (seqlens <= 258) & (seqlens >= 256)]

##### Remove chimeras - already done by core... but check to make sure
storage.mode(seqtab.filt) <- "integer"
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)
# Identified 0 bimeras out of 5722 input sequences.

#######################################################################
########### SAVE THE RDS FILE SO DON'T HAVE TO REPEAT ABOVE ###########
#######################################################################
saveRDS(seqtab.nochim, file="~/boxdrive/_Mrug/Diet/Microbiome_QC/raw/seqtab.nochim.rds")


##################### SKIP ABOVE, START HERE ############################
#########################################################################
### RELOAD THE SAVED INFO FROM HERE (if you have closed the project): ### 
#########################################################################
seqtab.nochim <- readRDS("./raw/seqtab.nochim.rds")


#########################################################################
################### ASSIGNING THE TAXONOMY WITH DADA2 ###################
#########################################################################

######################## SILVA Silva version 138 ########################
# Downloaded the Silva version 138 database here https://zenodo.org/record/3986799#.X6X7wFNKhBw
# import silva_nr99_v138_train_set.fa.gz
taxa <- dada2::assignTaxonomy(seqtab.nochim, "~/boxdrive/_Mrug/Diet/Microbiome_QC/reference_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
# Make species level assignments based on exact matching between ASVs and sequenced reference strains. 
# Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign 
# species to 16S gene fragments. Currently, species-assignment training fastas are available for the Silva 
# and RDP 16S databases. Download the silva_species_assignment_v138.fa.gz file.
#import silva_species_assignment_v138.fa.gz
taxa <- dada2::addSpecies(taxa, "~/boxdrive/_Mrug/Diet/Microbiome_QC/reference_db/silva_species_assignment_v138.fa.gz")

# Inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Export the file - this code removes the quotation marks 
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)

# Inspect the taxonomic assignments again:
taxon.print <- taxon # Removing sequence rownames for display only
rownames(taxon.print) <- NULL
head(taxon.print)

# Export the files
write.table(taxon,"~/boxdrive/_Mrug/Diet/R_analysis/_bin/silva_tax_table_NAs.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "~/boxdrive/_Mrug/Diet/R_analysis/_bin/dietw1_asv_table.txt",sep="\t",col.names=NA)

# Now export a separate set of files with the NAs removed - may want to use for some applications
# this replaces the NA with the nest highest classification

# Fix the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- paste("K", sep = "_", taxon$Kingdom[is.na(taxon$Phylum)])
taxon$Class[is.na(taxon$Class)] <- paste("P", sep = "_", taxon$Phylum[is.na(taxon$Class)])
taxon$Order[is.na(taxon$Order)] <- paste("C", sep = "_", taxon$Class[is.na(taxon$Order)])
taxon$Family[is.na(taxon$Family)] <- paste("O", sep = "_", taxon$Order[is.na(taxon$Family)])
taxon$Genus[is.na(taxon$Genus)] <- paste("F", sep = "_", taxon$Family[is.na(taxon$Genus)])
taxon$Species[is.na(taxon$Species)] <- paste("G", sep = "_", taxon$Genus[is.na(taxon$Species)])

# Inspect the taxonomic assignments again:
taxon.print <- taxon # Removing sequence rownames for display only
rownames(taxon.print) <- NULL
head(taxon.print)

# Export the file with NAs removed
write.table(taxon,"~/boxdrive/_Mrug/Diet/R_analysis/_bin/silva_tax_table_noNAs.txt",sep="\t",col.names=NA)
######################## SILVA complete ########################

######################## RDP trainset 18 ########################
# Downloaded the RDP trainset 18/release 11.5 here https://zenodo.org/record/4310151#.YtGgtuzMLX0
taxa <- dada2::assignTaxonomy(seqtab.nochim, "~/boxdrive/_Mrug/Diet/Microbiome_QC/reference_db/rdp_train_set_18.fa.gz", multithread=TRUE)
# Make species level assignments based on exact matching between ASVs and sequenced reference strains. 
taxa <- dada2::addSpecies(taxa, "~/boxdrive/_Mrug/Diet/Microbiome_QC/reference_db/rdp_species_assignment_18.fa.gz")

# Export the file - this code removes the quotation marks 
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)

# Export the files
write.table(taxon,"~/boxdrive/_Mrug/Diet/R_analysis/_bin/RDP_tax_table_NAs.txt",sep="\t",col.names=NA)

# Now export a separate set of files with the NAs removed
# Fix the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- paste("K", sep = "_", taxon$Kingdom[is.na(taxon$Phylum)])
taxon$Class[is.na(taxon$Class)] <- paste("P", sep = "_", taxon$Phylum[is.na(taxon$Class)])
taxon$Order[is.na(taxon$Order)] <- paste("C", sep = "_", taxon$Class[is.na(taxon$Order)])
taxon$Family[is.na(taxon$Family)] <- paste("O", sep = "_", taxon$Order[is.na(taxon$Family)])
taxon$Genus[is.na(taxon$Genus)] <- paste("F", sep = "_", taxon$Family[is.na(taxon$Genus)])
taxon$Species[is.na(taxon$Species)] <- paste("G", sep = "_", taxon$Genus[is.na(taxon$Species)])

# Inspect the taxonomic assignments again:
taxon.print <- taxon # Removing sequence rownames for display only
rownames(taxon.print) <- NULL
head(taxon.print)

# Export the file with NAs removed
write.table(taxon,"~/boxdrive/_Mrug/Diet/R_analysis/_bin/RDP_tax_table_noNAs.txt",sep="\t",col.names=NA)
######################## RDP complete ########################


###################################################################
################### CONSTRUCT PHYLOGENETIC TREE ################### 
###################################################################
# https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html
# https://github.com/benjjneb/dada2/issues/88

#Extract sequences from DADA2 output
sequences<-dada2::getSequences(seqtab.nochim)
names(sequences)<-sequences

# Run Sequence Alignment (MSA) using DECIPHER
alignment <- DECIPHER::AlignSeqs(DNAStringSet(sequences), anchor=NA)

# Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# Create distance matrix
dm <- dist.ml(phang.align)

# Perform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order

# Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)

# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

##### SAVE THE FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS
saveRDS(fitGTR, file="~/boxdrive/_Mrug/Diet/R_analysis//_bin/fitGTR.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
fitGTR <- readRDS("~/boxdrive/_Mrug/Diet/R_analysis/_bin/fitGTR.rds")


###################################################################
##################### CREATE PHYLOSEQ OBJECT ###################### 
###################################################################

###################### START HERE: IMPORT INTO PHYLOSEQ

### Load libraries ###
library(tidyverse)
library(phyloseq)
library(genefilter)
library(microbiome)
library(Biostrings)

### Set working directory ###
setwd("~/boxdrive/_Mrug/Diet") #desktop
# setwd("~/Box/_Mrug/Diet") #laptop

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("./R_analysis/_bin/dietw1_asv_table.txt",sep="\t",header=TRUE, row.names=1)
samples<-read.table("./R_analysis/_bin/mapping_file_SID.txt",sep="\t",header=T,row.names=1)
taxon<-read.table("./R_analysis/_bin/silva_tax_table_noNAs.txt",sep="\t",header=T,row.names=1)
fitGTR <- readRDS("./R_analysis/_bin/DietW1_fitGTR.rds")

# prepare phyloseq objects
OTU <- otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX <-tax_table(taxon)
TREE<-phy_tree(fitGTR$tree)
sampledata = sample_data(samples)

# combine
ps <- phyloseq(OTU, TAX, TREE, sampledata)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5722 taxa and 137 samples ]
# tax_table()   Taxonomy Table:    [ 5722 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5722 tips and 5720 internal nodes ]

rank_names(ps)
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family" 
# [6] "Genus"   "Species"

### Remove host contaminants, chloroplasts, mitochondria, and Eukaryota
ntaxa(ps) #5722
ps <- subset_taxa(ps, Kingdom !="Unknown")
ntaxa(ps) #5722; 0 removed
ps <- subset_taxa(ps, Family !="Mitochondria")
ntaxa(ps) #5719; 3 removed
ps <- subset_taxa(ps, Order !="Chloroplast")
ntaxa(ps) #5687; 32 removed

# check
sum(sample_sums(ps)) # 6789449 reads; 5687 taxa
get_taxa_unique(ps, "Kingdom") # okay
get_taxa_unique(ps, "Phylum") # okay
get_taxa_unique(ps, "Class") # okay
get_taxa_unique(ps, "Order") # okay
get_taxa_unique(ps, "Family") # okay
get_taxa_unique(ps, "Genus") # okay

# remove rare ASVs
# A -- The count value minmum threshold
# k -- The number of samples in which a taxa exceeded A
# A taxa is retained in the dataset if it exceeds the value A in at least k samples.
flist = genefilter::filterfun(kOverA(3, 5))
Diet <- filter_taxa(ps, flist, prune=TRUE)
ntaxa(Diet) #1254
(sum(sample_sums(ps)) - sum(sample_sums(Diet))) / sum(sample_sums(ps)) 
# 0.03295982 of the original sequences removed  

# mean, max and min of sample read counts
min(sample_sums(Diet)) #4087
mean(sample_sums(Diet)) #47924.6
max(sample_sums(Diet)) #91073

### add read count
reads_sample <- microbiome::readcount(Diet)
sample_data(Diet)$reads_sample <- reads_sample

# add a DNAStringSet object (to be accessed by the refseq function) to keep track of the ASV sequences
sequences <- Biostrings::DNAStringSet(taxa_names(Diet))
names(sequences) <- taxa_names(Diet)
Diet <- merge_phyloseq(Diet, sequences)
Diet
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5722 taxa and 137 samples ]
# tax_table()   Taxonomy Table:    [ 5722 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5722 tips and 5720 internal nodes ]
# refseq()      DNAStringSet:      [ 5722 reference sequences ]

### rename our ASVs to something more convenient in downstream analysis, while automatically retaining the corresponding unique sequence identifier
taxa_names(Diet) <- paste0("ASV", seq(ntaxa(Diet)))
view(taxa_names(Diet))

### save a key of ASV#'s and sequences
asv_key = as((refseq(Diet)), "data.frame")
colnames(asv_key )[1] <- "sequence"
rownames(asv_key) -> asv_key$ASV
write.table(asv_key,"./R_analysis/_bin/Diet_asv_key.txt",sep="\t",row.names=TRUE)

### now attach ASV#'s to the 4 taxon tables
# first read back in the tables because I didn't save as individual objects
Silva_noNAs<-read.table("./Microbiome_QC/silva_tax_table_noNAs.txt",sep="\t",header=T)
colnames(Silva_noNAs)[1] <- "sequence"
Silva<-read.table("./Microbiome_QC/silva_tax_table_NAs.txt",sep="\t",header=T)
colnames(Silva)[1] <- "sequence"
RDP_noNAs<-read.table("./Microbiome_QC/RDP_tax_table_noNAs.txt",sep="\t",header=T)
colnames(RDP_noNAs)[1] <- "sequence"
RDP<-read.table("./Microbiome_QC/RDP_tax_table_NAs.txt",sep="\t",header=T)
colnames(RDP)[1] <- "sequence"

### Now make a new tax table of filters ASVs
Silva_noNAs_tb <- dplyr::left_join(asv_key,Silva_noNAs) %>% select(.,-sequence) %>% column_to_rownames(.,"ASV")
Silva_tb <- dplyr::left_join(asv_key,Silva) %>% select(.,-sequence) %>% column_to_rownames(.,"ASV")
RDP_noNAs_tb <- dplyr::left_join(asv_key,RDP_noNAs) %>% select(.,-sequence) %>% column_to_rownames(.,"ASV")
RDP_tb <- dplyr::left_join(asv_key,RDP)%>% select(.,-sequence) %>% column_to_rownames(.,"ASV")


##################################################################
####### Export Cleaned OTU and Taxa Tables for Future Use ########
##################################################################

### export cleaned otu and taxa tables from phyloseq ### 
asv = as(otu_table(Diet), "matrix")
metadata = as(sample_data(Diet), "matrix") # need to delete the one low-read sample for use with phyloseq
tree = phy_tree(Diet)
write.table(asv,"./R_analysis/_bin/Diet_asv_table.txt",sep="\t",col.names=NA)
write.table(metadata,"./R_analysis/Microbiome_QC/Diet_map_seqreads.txt",sep="\t",col.names=NA)
ape::write.tree(tree, "./R_analysis/_bin/Diet.tree") #this tree will have the ASV# instead of the sequences

### export fasta refseq file
Diet %>%
  refseq() %>%
  Biostrings::writeXStringSet("./R_analysis/_bin/Diet_refset.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

### export taxon tables
# these taxon tables have the ASV# instead of the sequence and they also have the very rare taxa removed (see line 240)
write.table(Silva_noNAs_tb,"./R_analysis/_bin/Diet_Silva_noNAs_tax_table.txt",sep="\t",col.names=NA, row.names=TRUE)
write.table(Silva_tb,"./R_analysis/_bin/Diet_Silva_tax_table.txt",sep="\t",col.names=NA, row.names=TRUE)
write.table(RDP_noNAs_tb,"./R_analysis/_bin/Diet_RDP_noNAs_tax_table.txt",sep="\t",col.names=NA, row.names=TRUE)
write.table(RDP_tb,"./R_analysis/_bin/Diet_RDP_tax_table.txt",sep="\t",col.names=NA, row.names=TRUE)


##################################################################
####### Re-Import and Test theCleaned OTU and Taxa Tables ########
##################################################################

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("./R_analysis/_bin/Diet_asv_table.txt",sep="\t",header=TRUE, row.names=1)
#samples<-read.table("./Microbiome_QC/Diet_map_seqreads.txt",sep="\t",header=T,row.names=1)
samples<-read.table("./R_analysis/_bin/Diet_full_map.txt",sep="\t",header=T,row.names=1)
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
# sample_data() Sample Data:       [ 137 samples by 23 sample variables ]
# tax_table()   Taxonomy Table:    [ 1253 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1253 tips and 1251 internal nodes ]

# Remove the one sample with less than 5000 reads
ps <- subset_samples(ps, reads_sample>5000)

# mean, max and min of sample read counts
min(sample_sums(ps)) #24341
mean(sample_sums(ps)) #48246.93
max(sample_sums(ps)) #91073
sd(sample_sums(ps))/sqrt(136) #907.2953




