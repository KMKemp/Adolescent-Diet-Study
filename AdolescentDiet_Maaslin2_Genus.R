
# Set Working Directory ---------------------------------------------------

#setwd("~/boxdrive/_Mrug/Diet/") #desktop
setwd("~/Library/CloudStorage/Box-Box/_Mrug/Diet") #laptop

# Load Libraries ----------------------------------------------------------

library(Maaslin2)
library(phyloseq)
library('stringr')
library(tidyverse)
library(metagMisc)
library(compositions)
library(ggdendro)

# https://github.com/biobakery/biobakery/wiki/maaslin2


# Load data ---------------------------------------------------------------

# set files paths
otu <- read.table("./R_analysis/_bin/Diet_asv_table.txt",sep="\t",header=TRUE, row.names=1)
samples<-read.table("./R_analysis/_bin/Diet_full_map.txt",sep="\t",header=T,row.names=2)
taxa<-read.table("./R_analysis/_bin/Diet_Silva_tax_table.txt",sep="\t",header=T,row.names=1)
taxa <- as.matrix(taxa)

# make phyloseq object
Diet <- phyloseq(otu_table(t(otu), taxa_are_rows=TRUE), tax_table(taxa),
                 sampledata = sample_data(samples))


# Prep data and rarefy ----------------------------------------------------

# rarefy to minimum sample size
set.seed(5455)
rare_ps <- phyloseq::rarefy_even_depth(Diet)

metadata <- samples %>% dplyr::filter(reads_sample > 5000)
metadata <- metadata %>% mutate(ceth4cat = factor(ceth4cat, levels = c(1, 2, 3,4))) %>%
                         mutate(locale = factor(locale, levels = c(1, 2, 3)))
str(metadata)


# Summarize at the Phylum, Family, and Genus Levels -----------------------

# agglomerate at the Genes level
ps_genus <- tax_glom(rare_ps, taxrank = "Genus")
get_taxa_unique(ps_genus, "Genus") #225 genera
genus_df <- (as(otu_table(ps_genus), "matrix")) 
genus_df <- as.data.frame(genus_df) %>% rownames_to_column(.,"ASV")
genus_df_tax <- (as(tax_table(ps_genus), "matrix")) 
genus_df_tax <- as.data.frame(genus_df_tax) %>% rownames_to_column(.,"ASV")
genus_tbl <- dplyr::full_join(genus_df_tax,genus_df)
gen_ml2 <- genus_tbl %>% dplyr::select(Genus,9:144)
gen_ml2 <- as.data.frame(t(gen_ml2)) %>%
  `colnames<-`(.[1, ]) %>% .[-1, ] %>% mutate_if(is.character,as.numeric)


# Global Demographic Models -----------------------------------------------

#######################################################
######### Maaslin2- Demographics - Genus level ########
#######################################################

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/demographics", 
  fixed_effects = c("female", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  reference = c('ceth4cat,1;locale,1'),
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE) # already standardized
global_res <- fit_data[['results']]

# combine factor name and factor level
global_res$Variable <- paste(global_res$metadata, global_res$value)

# rename cols
global_res <- global_res %>%
  rename("Genus" = feature)


# add columns that tell us whether the p-value was less than 0.05
cov_coef_sig <- global_res %>%
  mutate(sig_q = ifelse(qval < 0.15, T, F), p_if_sig = ifelse(qval < 0.15, pval, NA), coef_if_sig = ifelse(qval < 0.15, coef, NA))

cov_coef_sig <- cov_coef_sig %>% dplyr::filter(!Genus %in% "Escherichia.Shigella") 

sig_taxa <-cov_coef_sig %>%
  group_by(Genus) %>%
  filter(meanRow(p_if_sig) > 0 ) %>%
  pull(Genus)


sig_cov <- cov_coef_sig%>%
  group_by(Variable) %>%
  filter(meanRow(p_if_sig) > 0 ) %>%
  pull(Variable)

cov_coef_sig_f <- cov_coef_sig %>% dplyr::filter(Genus %in% sig_taxa) 
cov_coef_sig_f <- cov_coef_sig_f %>% dplyr::filter(Variable %in% sig_cov ) 

cov_coef_sig_f$Genus<-gsub(".", " ", cov_coef_sig_f$Genus, fixed=TRUE)
cov_coef_sig_f$Genus<-gsub("R 7", "R-7", cov_coef_sig_f$Genus, fixed=TRUE)
cov_coef_sig_f$Genus<-gsub("UCG ", "UCG-", cov_coef_sig_f$Genus, fixed=TRUE)
cov_coef_sig_f$Genus<-gsub("X ", "[", cov_coef_sig_f$Genus, fixed=TRUE)
cov_coef_sig_f$Genus<-gsub("  ", "] ", cov_coef_sig_f$Genus, fixed=TRUE)

cov_coef_sig_f2 <-cov_coef_sig_f %>% dplyr::filter(sig_q == TRUE)

genus_tbl <- genus_tbl %>% select(Phylum, Class, Order, Family, Genus)
global_gen_tbl <- dplyr::right_join(genus_tbl, cov_coef_sig_f2)
global_gen_tbl  <-arrange(global_gen_tbl  , Variable, Phylum, Class, Order, Family)


write.csv(global_gen_tbl,"./R_analysis/Maaslin2/Genus_level/demographics/Maaslin2_NEGBIN_genus_global_FDR.15_10.6.2022.csv") 



# REAP Models -------------------------------------------------------------

#######################################################
############ Maaslin2- REAP - Genus level #############
#######################################################

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reapx", 
  fixed_effects = c("c1reapx", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reapx_res <- fit_data[['results']]
c1reapx_res <- c1reapx_res %>% dplyr::filter(metadata == "c1reapx")


fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap1", 
  fixed_effects = c("c1reap1", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap1_res <- fit_data[['results']]
c1reap1_res <- c1reap1_res %>% dplyr::filter(metadata == "c1reap1")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap2", 
  fixed_effects = c("c1reap2", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap2_res <- fit_data[['results']]
c1reap2_res <- c1reap2_res %>% dplyr::filter(metadata == "c1reap2")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap3", 
  fixed_effects = c("c1reap3", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap3_res <- fit_data[['results']]
c1reap3_res <- c1reap3_res %>% dplyr::filter(metadata == "c1reap3")


fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap4", 
  fixed_effects = c("c1reap4", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap4_res <- fit_data[['results']]
c1reap4_res <- c1reap4_res %>% dplyr::filter(metadata == "c1reap4")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap4a", 
  fixed_effects = c("c1reap4a", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap4a_res <- fit_data[['results']]
c1reap4a_res <- c1reap4a_res %>% dplyr::filter(metadata == "c1reap4a")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap5", 
  fixed_effects = c("c1reap5", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap5_res <- fit_data[['results']]
c1reap5_res <- c1reap5_res %>% dplyr::filter(metadata == "c1reap5")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap6", 
  fixed_effects = c("c1reap6", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap6_res <- fit_data[['results']]
c1reap6_res <- c1reap6_res %>% dplyr::filter(metadata == "c1reap6")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap7", 
  fixed_effects = c("c1reap7", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap7_res <- fit_data[['results']]
c1reap7_res <- c1reap7_res %>% dplyr::filter(metadata == "c1reap7")


fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap8", 
  fixed_effects = c("c1reap8", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap8_res <- fit_data[['results']]
c1reap8_res <- c1reap8_res %>% dplyr::filter(metadata == "c1reap8")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap9", 
  fixed_effects = c("c1reap9", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap9_res <- fit_data[['results']]
c1reap9_res <- c1reap9_res %>% dplyr::filter(metadata == "c1reap9")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap10", 
  fixed_effects = c("c1reap10", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap10_res <- fit_data[['results']]
c1reap10_res <- c1reap10_res %>% dplyr::filter(metadata == "c1reap10")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap11", 
  fixed_effects = c("c1reap11", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap11_res <- fit_data[['results']]
c1reap11_res <- c1reap11_res %>% dplyr::filter(metadata == "c1reap11")

fit_data = Maaslin2(
  input_data = gen_ml2, 
  input_metadata = metadata, 
  output = "./R_analysis/Maaslin2/Genus_level/genusc1reap12", 
  fixed_effects = c("c1reap12", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
c1reap12_res <- fit_data[['results']]
c1reap12_res <- c1reap12_res %>% dplyr::filter(metadata == "c1reap12")

#######################################################
###### combine result tables for all reap items #######
#######################################################

reap_coef <- rbind(c1reapx_res, c1reap1_res,c1reap2_res,c1reap3_res,c1reap4_res,c1reap4a_res,
                     c1reap5_res, c1reap6_res, c1reap7_res, c1reap8_res, c1reap9_res,
                     c1reap10_res, c1reap11_res, c1reap12_res)
# rename cols
reap_coef <- reap_coef %>%
  rename("Genus" = feature, "Variable"= metadata)

# add columns that tell us whether the p-value was less than 0.05
reap_coef_sig <- reap_coef %>%
  mutate(sig_q = ifelse(qval < .15, T, F), p_if_sig = ifelse(qval < .15, pval, NA), r_if_sig = ifelse(qval < .15, coef, NA))


sig_taxa <- reap_coef_sig %>%
  group_by(Genus) %>%
  filter(!is.na(p_if_sig)) %>%
  pull(Genus)

sig_reaps <- reap_coef_sig %>%
  group_by(Variable) %>%
  filter(!is.na(p_if_sig)) %>%
  #filter(str_detect(Variable, "c1reap")) %>%
  pull(Variable)


reap_coef_sig_f <- reap_coef_sig %>% dplyr::filter(Genus %in% sig_taxa) 
reap_coef_sig_f <- reap_coef_sig_f %>% dplyr::filter(Variable %in% sig_reaps) 

reap_coef_sig_f$Genus<-gsub(".", " ", reap_coef_sig_f$Genus, fixed=TRUE)
reap_coef_sig_f$Genus<-gsub("R 7", "R-7", reap_coef_sig_f$Genus, fixed=TRUE)

reap_coef_sig_f2 <- reap_coef_sig_f %>% filter(sig_q == TRUE)

genus_tbl <- genus_tbl %>% select(Phylum, Class, Order, Family, Genus)
reap_gen_tbl <- dplyr::right_join(genus_tbl, reap_coef_sig_f2)
reap_gen_tbl  <-arrange(reap_gen_tbl , Variable, Phylum, Class, Order, Family)

write.csv(reap_gen_tbl,"./R_analysis/Maaslin2/Genus_level/maaslin2_genus_reap_FDR.15_10.5.2022.csv") 


# Demographics Heatmap ----------------------------------------------------

test <- cov_coef_sig_f %>% select(Genus, Variable, coef) %>% 
  pivot_wider(names_from = Genus, values_from = coef) %>% column_to_rownames("Variable")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)

demo_order <- c("Phylum","c1zbmi c1zbmi" ,"SES SES" ,"female female" ,  "ceth4cat 2" ,  "ceth4cat 4" ,
"ceth4cat 3" ,  "locale 2"  ,    "locale 3")  

test <- cov_coef_sig_f %>% select(Genus, Variable, coef) %>% 
  pivot_wider(names_from = Variable, values_from = coef) %>% column_to_rownames("Genus")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)

list <- print(order.dendrogram(dendro)) 
Taxa.order <- print(row.names(test[c(list),]))
Taxa.order


demo_hm <- cov_coef_sig_f   %>% 
  ggplot(aes(Variable, Genus, fill=coef, label=round(coef_if_sig,2))) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Masslin2\nCorrelation", title="", subtitle="") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0072B2",high="#D55E00", limits=c(-4,4)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0), limits = demo_order,
                   labels=c("female female"="female",  "c1zbmi c1zbmi"="zBMI", "SES SES"="SES" ,
                            "locale 2" = "suburb","locale 3" = "city","ceth4cat 2" = "Black",
                            "ceth4cat 3" = "Hispanic","ceth4cat 4" = "other\nminority")) +
  scale_y_discrete(expand=c(0,0),
                   limits = Taxa.order) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95)) +
  theme(axis.text = element_text(size = 11, color = "black"))

demo_hm

# export eps as 1250 x 450

# Reap/Diet Heatmap -------------------------------------------------------

test <- reap_coef_sig_f %>% select(Genus, Variable, coef) %>% 
  pivot_wider(names_from = Genus, values_from = coef) %>% column_to_rownames("Variable")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro)) 
reap.order <- print(row.names(test[c(list),]))
reap.order <- c("Phylum","Family","c1reapx" ,  "c1reap4" , "c1reap4a", "c1reap5" ,"c1reap11", "c1reap12", "c1reap9" ,"c1reap8" , 
"c1reap2" , "c1reap6",  "c1reap7" , "c1reap1" , "c1reap3" )


test <- reap_coef_sig_f %>% select(Genus, Variable, coef) %>% 
  pivot_wider(names_from = Variable, values_from = coef) %>% column_to_rownames("Genus")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro))
list <- print(order.dendrogram(dendro)) 
taxa.order <- print(row.names(test[c(list),]))
taxa.order




reap_hm<- reap_coef_sig_f   %>% 
  ggplot(aes(Variable, Genus, fill=coef, label=round(r_if_sig,2))) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Masslin2\nCorrelation", title="", subtitle="") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0072B2",high="#D55E00", limits=c(-2.0,2.0)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0), limits = reap.order,
                   labels=c("c1reap1"="skipping\nbreakfast", "c1reap2"="eating out", "c1reap3"="grains",
                            "c1reap4"="fruits" , "c1reap4a" = "vegetables","c1reap5" = "dairy/cheese",
                            "c1reap6" = "meat",  "c1reap7" = "processed\nmeat","c1reap8" = "fried foods",
                            "c1reap9" = "fatty snacks", "c1reap11" = "sweets",
                            "c1reap12" = "soads\nsugary drinks","c1reapx" = "overall\ndiet quality"))+
  scale_y_discrete(expand=c(0,0),
                   limits = taxa.order) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95)) +
  theme(axis.text = element_text(size = 11, color = "black"))
reap_hm




plot_grid(demo_hm, reap_hm,
          #provide labels for plots
          labels = c('A', 'B'), 
          #ensure plots are aligned on the horizontal axis
          align = "hv",
          axis = c("lr"),
          ncol = 1, nrow = 2,
          greedy = TRUE,
          rel_heights = c(3,2))
