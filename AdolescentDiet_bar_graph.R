

### Set working directory ###
setwd("~/boxdrive/_Mrug/Diet") #desktop
#setwd("~/Box/_Mrug/Diet") #laptop

### Load libraries ###
library(tidyverse)
library(phyloseq)
library(ggplot2)


####################################################################
################### Import into Phyloseq for QC ####################
####################################################################

# set files paths
otu <- read.table("./R_analysis/_bin/Diet_asv_table.txt",sep="\t",header=TRUE, row.names=1)
samples<-read.table("./R_analysis/_bin/Diet_full_map.txt",sep="\t",header=T,row.names=2)
taxa<-read.table("./R_analysis/_bin/Diet_Silva_tax_table.txt",sep="\t",header=T,row.names=1)
taxa <- as.matrix(taxa)

# make phyloseq object
Diet <- phyloseq(otu_table(t(otu), taxa_are_rows=TRUE), tax_table(taxa),
                 sampledata = sample_data(samples))

# Prep data and rarefy ----------------------------------------------------

# filter out ASVs that are found in less than 10% of samples
#Diet <- metagMisc::phyloseq_filter_prevalence(Diet, prev.trh = 0.25, abund.trh = NULL)

# rarefy to minimum sample size
set.seed(5455)
rare_ps <- phyloseq::rarefy_even_depth(Diet)


# Transform to relative abundance
Diet_norm <- transform_sample_counts(rare_ps, function(x) 100 * x/sum(x))

# Prepare to make a barchart at order level
Diet_bc <- Diet_norm  %>%
  tax_glom(taxrank = "Phylum") %>% # agglomerate taxa at order level
  psmelt() %>%                    # Melt phyloseq object to long format for producing graphics with ggplot2
  filter(Abundance >4.0)  %>%    # Filter out orders below 1% in each sample
  arrange(desc(Phylum))

# Check how many Orders 
unique(Diet_bc$Phylum) #5

# Sum remaining taxa with a relative abundance < 1% and make a new dataframe
Remainders <- (Diet_bc) %>%
  dplyr::group_by(SID) %>% 
  dplyr::summarise(Abundance = (100-sum(Abundance))) %>% 
  as.data.frame()
Remainders$Phylum<-"Phyla < 4% RA"


# Compile dataframes
Diet_bc <- full_join(Diet_bc,Remainders)
Diet_bc$Phylum <- as.factor(Diet_bc$Phylum)

# lock in Order level order
Diet_bc$Phylum  <- factor(Diet_bc$Phylum, 
                             levels = c("Phyla < 4% RA","Verrucomicrobiota","Proteobacteria","Actinobacteriota","Bacteroidota",
                                       "Firmicutes"))

arrange(Diet_bc,Abundance, by_group = "Phylum")

Diet_bc <- (Diet_bc %>% arrange(desc(Abundance)) %>% arrange(desc(Phylum))) 
Diet_bc$SID <- fct_inorder(as.factor(Diet_bc$SID)) 

write.csv(Diet_bc , file = "Diet_bc_phylum.csv")

# Function for making a custom color palette for plotting
library(RColorBrewer)
# Function for plotting colors side-by-side from here: https://cran.r-project.org/web/packages/colorspace/index.html
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

colorBlindGrey8   <- c("#999999",  "#7e0018",
                       "#db6d00", "#00489e", "#005745", "#450270")
pal(colorBlindGrey8)




# Plot barchart for Orders with samples labeled
p1<- ggplot(Diet_bc, aes(x = as.factor(SID), y = Abundance, fill = Phylum))+
  geom_bar(stat = "identity", colour = "black",size=0.5,width=1)+
  scale_fill_manual(values = colorBlindGrey8 )+
  labs(x = " ",y = "Phylum Relative Abundance (% RA)")+
  theme(plot.title = element_text(hjust = 0.5, size=14,face="bold"))+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  theme(axis.text.y = element_text(size=11), axis.title.y=element_text(size=12))+
  theme(legend.text=element_text(size=11), legend.title=element_text(size=12,face="bold"))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=0.75, linetype="solid"))+
  theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 10))+
  scale_x_discrete(expand = c(0,0))
p1



