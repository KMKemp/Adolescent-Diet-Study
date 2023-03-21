
# Set Working Directory ---------------------------------------------------

setwd("~/boxdrive/_Mrug/Diet/") #desktop
#setwd("~/Box/_Mrug/Diet") #laptop


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


# Summarize at the Phylum Level -----------------------------------------

# agglomerate at the Phylym level
ps_phy <- tax_glom(rare_ps, taxrank = "Phylum", NArm = TRUE)
get_taxa_unique(ps_phy , "Phylum") #65 families
df_phy <- (as(otu_table(ps_phy), "matrix"))  
df_phy <- as.data.frame(df_phy) %>% rownames_to_column(.,"ASV")
df_phy_tax <- (as(tax_table(ps_phy), "matrix"))  
df_phy_tax <- as.data.frame(df_phy_tax) %>% rownames_to_column(.,"ASV")
phy_tbl <- dplyr::full_join(df_phy_tax,df_phy)
phy_ml2 <- phy_tbl %>% dplyr::select(Phylum,9:144)
phy_ml2 <- as.data.frame(t(phy_ml2)) %>%
  `colnames<-`(.[1, ]) %>% .[-1, ] %>% mutate_if(is.character,as.numeric)



# Global Demographic Model ------------------------------------------------

#######################################################
######## Maaslin2- Demographics - Phylum level ########
#######################################################

fit_data = Maaslin2(
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/demographic", 
  fixed_effects = c("female", "female", "c1zbmi", "ceth4cat","locale","SES"), 
  reference = c('ceth4cat,1;locale,1'),
  analysis_method = "NEGBIN", 
  transform = "NONE",
  normalization = "TMM",
  min_prevalence = 0.75,
  min_abundance =0,
  max_significance = 0.15,
  standardize = FALSE, # already standardized
  plot_heatmap = FALSE,
  plot_scatter = FALSE) 
global_res <- fit_data[['results']]

# combine factor name and factor level
global_res$Variable <- paste(global_res$metadata, global_res$value)

# rename cols
global_res <- global_res %>%
  rename("Phylum" = feature)

# add columns that tell us whether the p-value was less than 0.05
cov_coef_sig <- global_res %>%
  mutate(sig_q = ifelse(qval < 0.15, T, F), p_if_sig = ifelse(qval < 0.15, pval, NA), r_if_sig = ifelse(qval < 0.15, coef, NA))

sig_taxa <-cov_coef_sig %>%
  group_by(Phylum) %>%
  filter(!is.na(p_if_sig)) %>%
  pull(Phylum)

sig_cov <- cov_coef_sig%>%
  group_by(Variable) %>%
  filter(!is.na(p_if_sig)) %>%
  pull(Variable)

#cov_coef_sig_f <- cov_coef_sig %>% dplyr::filter(Phylum %in% sig_taxa) 
cov_coef_sig_f <- cov_coef_sig %>% dplyr::filter(Variable %in% sig_cov ) 

cov_coef_sig_f2 <-cov_coef_sig_f %>% dplyr::filter(sig_q == TRUE)

phy_tbl <- phy_tbl %>% select(Phylum, Class, Order, Family, Genus)
global_phy_tbl <- dplyr::right_join(phy_tbl, cov_coef_sig_f2)
global_phy_tbl  <-arrange(global_phy_tbl  , Variable, Phylum, Class, Order, Family)


write.csv(global_phy_tbl,"~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/demographic/Maaslin2_NEGBIN_Phy_global_FDR.15_9.20.2022.csv") 


# REAP Items at the Phylum Level ------------------------------------------

#######################################################
############### Maaslin2 - Phylum level ###############
############## results for reap variables #############
#######################################################

fit_data = Maaslin2(
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_diet_quality", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_skipping_breakfast", 
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
  input_data = phy_ml2,
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_eatout", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_wholegrains", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_fruit", 
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
  input_data = phy_ml2,
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_vegetable", 
  fixed_effects = c("c1reap4a","female", "c1zbmi", "ceth4cat","locale","SES"), 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_dairy", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_meat", 
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
  input_data = phy_ml2,
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_processedmeat", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_friedfoods", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_fattysnacks", 
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
  input_data = phy_ml2,
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_fats_oils", 
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
  input_data = phy_ml2, 
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_sweets", 
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
  input_data = phy_ml2,
  input_metadata = metadata, 
  output = "~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/phylum_sodas_sugarydrinks", 
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
  rename("Phylum" = feature, "Variable"= metadata)

# add columns that tell us whether the p-value was less than 0.05
reap_coef_sig <- reap_coef %>%
  mutate(sig_q = ifelse(qval < .15, T, F), p_if_sig = ifelse(qval < .15, pval, NA), r_if_sig = ifelse(qval < .15, coef, NA))

sig_taxa <- reap_coef_sig %>%
  group_by(Phylum) %>%
  filter(!is.na(p_if_sig)) %>%
  pull(Phylum)

sig_cov <- reap_coef_sig %>%
  group_by(Variable) %>%
  filter(!is.na(p_if_sig)) %>% 
  pull(Variable)

#reap_coef_sig_f <- reap_coef_sig %>% dplyr::filter(Phylum%in% sig_taxa) 
reap_coef_sig_f <- reap_coef_sig %>% dplyr::filter(Variable %in% sig_cov) 

reap_coef_sig_f2 <- reap_coef_sig_f %>% filter(sig_q == TRUE)

phy_tbl <- phy_tbl %>% select(Phylum, Class, Order, Family, Genus)
reap_phy_tbl <- dplyr::right_join(phy_tbl, reap_coef_sig_f2)
reap_phy_tbl  <-arrange(reap_phy_tbl, Variable, Phylum, Class, Order, Family)

write.csv(reap_phy_tbl,"~/boxdrive/_Mrug/Diet/R_analysis/Maaslin2/Phylum_level/Maaslin2_NEGBIN_Phylum_reap_FDR.15_10.7.2022.csv") 


# Demographics Heatmap ----------------------------------------------------

test <- cov_coef_sig_f %>% select(Phylum, Variable, coef) %>% 
  pivot_wider(names_from = Phylum, values_from = coef) %>% column_to_rownames("Variable")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro))
demo.order <- print(row.names(test[c(list),]))

demo.order <- c( "locale 3"   ,  "locale 2" ,   "ceth4cat 3",    "ceth4cat 4" ,   "ceth4cat 2" ,     "SES SES" ,  "c1zbmi c1zbmi"   )
demo.order <- rev(demo.order)


test <- cov_coef_sig_f %>% select(Phylum, Variable, coef) %>% 
  pivot_wider(names_from = Variable, values_from = coef) %>% column_to_rownames("Phylum")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro))
taxa.order <- print(row.names(test[c(list),]))



p_demo_hm <- cov_coef_sig_f   %>% 
  ggplot(aes(Variable, Phylum, fill=coef, label=round(r_if_sig,2))) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Masslin2\nCorrelation", title="", subtitle="") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0072B2",high="#D55E00", limits=c(-2.5,2.5)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0), limits = demo.order,
                   labels=c("female female"="female",  "c1zbmi c1zbmi"="zBMI", "SES SES"="SES" ,
                            "locale 2" = "suburb","locale 3" = "city","ceth4cat 2" = "Black",
                            "ceth4cat 3" = "Hispanic","ceth4cat 4" = "other\nminority")) +
  scale_y_discrete(expand=c(0,0),
                   limits = taxa.order) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95)) +
  theme(axis.text = element_text(size = 11, color = "black"))
p_demo_hm

# export eps as 830 x 430

# Reap/Diet Heatmap -------------------------------------------------------

test <- reap_coef_sig_f %>% select(Phylum, Variable, coef) %>% 
  pivot_wider(names_from = Phylum, values_from = coef) %>% column_to_rownames("Variable")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro))
reap.order <- print(row.names(test[c(list),]))


reap.order <- c(   "c1reap7" , "c1reap1", "c1reap2" ,"c1reap4a", "c1reap4") 
reap.order <- rev(reap.order)

test <- reap_coef_sig_f %>% select(Phylum, Variable, coef) %>% 
  pivot_wider(names_from = Variable, values_from = coef) %>% column_to_rownames("Phylum")
test <- as.matrix(test)
dendro <- as.dendrogram(hclust(d = dist(x = test)))
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro_plot)
list <- print(order.dendrogram(dendro))
taxa.order <- print(row.names(test[c(list),]))




p_reap_hm<- reap_coef_sig_f   %>% 
  ggplot(aes(Variable, Phylum, fill=coef, label=round(r_if_sig,2))) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Masslin2\nCorrelation", title="", subtitle="") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0072B2",high="#D55E00", limits=c(-1,1)) +
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
p_reap_hm

plot_grid(p_demo_hm, p_reap_hm, demo_hm, reap_hm, #the "demo_hm" and "reap_hm" plots are from the family level analyses... need to run that file first to create them
          #provide labels for plots
          labels = c('A', 'B', 'C', 'D'), 
          #ensure plots are aligned on the horizontal axis
          align = "vh",
          axis = c("lb"),
          ncol = 2, nrow = 2,
          greedy = FALSE,
          rel_widths = c(1, 1,4,4),
          rel_heights = c(1, 1,5,4.8))
