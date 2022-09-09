#libraries to load
library(qiime2R)
library(tidyverse)
library("stringr")  
library(ggplot2)
library("dplyr")
library(ggpubr)
library(rstatix)
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
#library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
#library(data.table) # alternative to data.frame
library(dplyr) # data handling  
library(picante)
library(phyloseq)

#load metadata file
metadata = read.table("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/metadata.txt", sep='\t', header=TRUE)

#read and parse taxonomy data
taxonomy<-read_qza('/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/taxonomy.qza')
head(taxonomy$data)
#parse taxonomy data into different levels
taxonomy_parsed<-parse_taxonomy(taxonomy$data)
head(taxonomy_parsed)

#read in all of the alpha diversity metrics that Koh lab uses 
#alpha-diversity: pielou evenness, faith, observed_features, shannon, 
shannon<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/shannon_vector.qza")
evenness<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/evenness_vector.qza")
observed_features<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/observed_features_vector.qza")
faith<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/faith_pd_vector.qza")

shannon_data = shannon$data
evenness_data = evenness$data
observed_features_data = observed_features$data
faith_data = faith$data

df_alpha <- data.frame (samples = rownames(shannon_data), Antibiotics = metadata$group, shannon  = shannon_data$shannon_entropy,
                        evenness = evenness_data$pielou_evenness, observed_feature = observed_features_data$observed_features, faith_pd = faith_data$faith_pd)

#create a phyloseq object
#requirements are the table.qza, rooted-tree.qza, taxonomy.qza and optional (metadata file)
physeq<-qza_to_phyloseq(
  features="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/table.qza",
  tree="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/rooted-tree.qza",
  taxonomy="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/taxonomy.qza",
  metadata = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/metadata.txt"
)

print(physeq)
#large difference in the number of reads (below)
summary(sample_sums(physeq))

#plot the rarifaction curve below to check how has the richness captured in the sequencing effort
otu_tab <- t(abundances(physeq))
p <- vegan::rarecurve(otu_tab, 
                      step = 50, label = FALSE, 
                      sample = min(rowSums(otu_tab), 
                                   col = "blue", cex = 0.6))
p
#y = number of species
#x = sample size
#Note: quality control script did not change phyloseq data
#this must be already incorporated into Jiwoong's pipeline

#biases in alpha diversity statistics
#https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full

#normalize data (spread in sample reads)
set.seed(9242)  # This will help in reproducing the filtering and nomalisation. 
#sample.size is based on the lowest amount of reads
#rarification randomly removes reads so that the read depths are the same
ps0.rar <- rarefy_even_depth(physeq, sample.size = 119)
#plot to see phylums present in samples after rarification
p.rar <- plot_taxa_prevalence(ps0.rar, "Phylum")
p.rar
#plot to see phylums present in samples with raw data
p.raw <- plot_taxa_prevalence(physeq, "Phylum")
p.raw

#faith pd code for rarified data
#get otu info
ps0.rar.asvtab <- as.data.frame(ps0.rar@otu_table)
#get tree info
ps0.rar.tree <- ps0.rar@phy_tree
# check whether tree is rooted or not. should be rooted if using above code!
ps0.rar@phy_tree
norm_pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).

#get alpha diversities for rarified data
norm_shannon <- alpha(ps0.rar, index = "diversity_shannon")
norm_evenness <- alpha(ps0.rar, index = "evenness_pielou")
norm_observed_feature<- alpha(ps0.rar, index = "observed")
#create table with the same labels as for the unnormalized data
norm_alpha_df <- data.frame (Antibiotics = df_alpha$Antibiotics, shannon = norm_shannon$diversity_shannon, evenness = norm_evenness$evenness_pielou, 
                             observed_feature = norm_observed_feature$observed, faith_pd = norm_pd$PD
)

#KW calculation for box plot (added to the figure)
KW_pval <- function(df, alpha_method) {
  if (alpha_method == "shannon") {
    pwc2 <- df %>% wilcox_test(shannon ~ Antibiotics, p.adjust.method = "bonferroni")
    res.kruskal <- df %>% kruskal_test(shannon ~ Antibiotics)
  }
  else if (alpha_method == "evenness") {
    pwc2 <- df %>% wilcox_test(evenness ~ Antibiotics, p.adjust.method = "bonferroni")
  } else if (alpha_method == "faith_pd") {
    pwc2 <- df %>% wilcox_test(faith_pd ~ Antibiotics, p.adjust.method = "bonferroni")
  } 
  else {
    pwc2 <- df %>% wilcox_test(observed_feature ~ Antibiotics, p.adjust.method = "bonferroni")
  }
  return(pwc2)
}

#KW calculation for box plot (added to the figure)
KW_all_test <- function(df, alpha_method) {
  if (alpha_method == "shannon") {
    res.kruskal <- df %>% kruskal_test(shannon ~ Antibiotics)
  }
  else if (alpha_method == "evenness") {
    res.kruskal <- df %>% kruskal_test(evenness ~ Antibiotics)
  } else if (alpha_method == "faith_pd") {
    res.kruskal <- df %>% kruskal_test(faith_pd ~ Antibiotics)
  } else {
    res.kruskal <- df %>% kruskal_test(observed_feature ~ Antibiotics)
  }
  return(res.kruskal)
}

#boxplot with KW Calculations for alpha diversity
boxplot_alpha <- function(df, alpha_method, plot_name) { 
  par(mar=c(.05,.000000000000000000001,.05,.05))
  #pwc2 <- df %>% 
  #  wilcox_test(shannon ~ Antibiotics, p.adjust.method = "bonferroni")
  pwc2 <- KW_pval(df, alpha_method)
  res.kruskal <- KW_all_test(df, alpha_method)
  pwc2 <- pwc2 %>% add_xy_position(x = "Antibiotics")
  p <- ggboxplot(df_alpha, x = "Antibiotics", y = alpha_method, fill = "Antibiotics", palette = 'jco') +
    stat_pvalue_manual(pwc2, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(res.kruskal, detailed = TRUE),
      caption = get_pwc_label(pwc2))
  p <- p + rotate_x_text()
  print(p)
  ggsave(path = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/alpha/", filename =  plot_name, device='tiff', dpi=300)
  return(p)
}

#normalized
boxplot_alpha( df_alpha, "shannon", "unnormalized_shannon_alpha")
boxplot_alpha( df_alpha, "faith_pd", "unnormalized_faith_alpha")
boxplot_alpha( df_alpha, "evenness", "unnormalized_evenness_alpha")
boxplot_alpha( df_alpha, "observed_feature", "unnormalized_observed_features_alpha")
#unnormalized
boxplot_alpha(norm_alpha_df , "shannon", "normalized_shannon_alpha")
boxplot_alpha( norm_alpha_df , "faith_pd", "normalized_faith_alpha")
boxplot_alpha( norm_alpha_df , "evenness", "normalized_evenness_alpha")
boxplot_alpha( norm_alpha_df , "observed_feature", "normalized_observed_features_alpha")

#read all beta diversity metrics that Koh lab uses
#beta diversity: BC, jaccard, unweighted Unifrac, weighted Unifrac
#read in all of the beta diversity metrics that Koh lab uses 
BC<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/bray_curtis_distance_matrix.qza")
BC_pcoa<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/bray_curtis_pcoa_results.qza")
jaccard<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/jaccard_distance_matrix.qza")
jaccard_pcoa<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/jaccard_pcoa_results.qza")
uw_un<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/unweighted_unifrac_distance_matrix.qza")
uw_un_pcoa<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/unweighted_unifrac_pcoa_results.qza")
w_un<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/weighted_unifrac_distance_matrix.qza")
w_un_pcoa<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/weighted_unifrac_pcoa_results.qza")

#plot pcoas for output files 
#Note: for the pcoa files, the eigenvals have already been calculated, so just use a scatterplot
#this function is for the previously calcualted pcoa files
#requires the loaded metadata file for colors
pcoa_beta <- function(df_pcoa, metadata_group, plot_name) { 
  df_pcoa = df_pcoa$data
  df_pcoa$Eigvals
  labx = paste("PCA1", as.character(df_pcoa$ProportionExplained[1,1]))
  laby = paste("PCA2", as.character(df_pcoa$ProportionExplained[1,2]))
  matrix <- df_pcoa$Vectors
  matrix["Treatment"] <- metadata_group
  p <- ggplot(matrix, aes(x=PC1, y=PC2, color= Treatment)) + geom_point() + labs(title = "PCoA", x = labx, y=laby)
  ggsave(path = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/beta/", filename =  plot_name, device='tiff', dpi=300)
  return(p)
  }

pcoa_beta(BC_pcoa, metadata$group, "precalculated_pcoa_BC")
pcoa_beta(jaccard_pcoa, metadata$group, "precalculated_pcoa_jaccard")
pcoa_beta(uw_un_pcoa, metadata$group, "precalculated_pcoa_uw_unifrac")
pcoa_beta(w_un_pcoa, metadata$group, "precalculated_pcoa_w_unifrac")

#for rarified data. Note: this can also be used for non rarfied data
# if we remove OTUs that are detected atleast 10 times in 5% of the samples
## we canreduce the sparsity considerably. 
#ps0.rar.filtered <- core(ps0.rar, detection = 10, prevalence = 0.05)
#summarize_phyloseq(ps0.rar.filtered)
#run quality control script - created by XZ
ps0.rar.qc <- qc(ps0.rar)

pcoa_beta_phylo_object <- function(qc_df, distance, metadata_group, plot_name) { 
  ps_ob_o = phyloseq::ordinate(qc_df, method="MDS", distance= distance)
  d2 = data.frame(x=ps_ob_o$vectors[,1], y=ps_ob_o$vectors[,2], group=metadata_group)
  p <- ggplot(d2, aes(x=x, y=y, color=group))+ geom_point() +labs(title = "PCoA", x= "PCA1", y="PCA2")
  ggsave(path = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/beta/", filename =  plot_name, device='tiff', dpi=300)
  return(p+stat_ellipse())
}

pcoa_beta_phylo_object(ps0.rar.qc, "unifrac", metadata$group, "rarified_un_unifrac")
pcoa_beta_phylo_object(ps0.rar.qc, "wunifrac", metadata$group, "rarified_w_unifrac")
pcoa_beta_phylo_object(ps0.rar.qc, "jaccard", metadata$group, "rarified_jaccard")
pcoa_beta_phylo_object(ps0.rar.qc, "bray", metadata$group, "rarified_braycurtis")




