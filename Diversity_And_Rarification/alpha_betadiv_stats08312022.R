#helpful links
#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#https://www.r-bloggers.com/2021/08/how-to-perform-tukey-hsd-test-in-r/
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/


#packages to install
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
install.packages("stringr")
install.packages("dplyr")    
install.packages("ggplot2") 

#libraries to load
library(qiime2R)
library(tidyverse)
library("stringr")  
library(ggplot2)
library("dplyr")
library(ggpubr)
library(rstatix)

#load metadata file
metadata = read.table("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/metadata.txt", sep='\t', header=TRUE)

#read and parse taxonomy data
taxonomy<-read_qza('/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/taxonomy.qza')
head(taxonomy$data)
#parse taxonomy data into different levels
taxonomy_parsed<-parse_taxonomy(taxonomy$data)
head(taxonomy_parsed)

#create a phyloseq object
#requirements are the table.qza, rooted-tree.qza, taxonomy.qza and optional (metadata file)
physeq<-qza_to_phyloseq(
  features="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/table.qza",
  tree="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/rooted-tree.qza",
  taxonomy="/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/taxonomy.qza",
  metadata = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/metadata.txt"
)

#Shannon diversity: how difficult it is to predict the identity of a randomly chosen individual
#read in shannon diversity
shannon<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/shannon_vector.qza")

#venn diagram control 
#number of shannon results should match the metadata number
#if not, check this again
shannon<-shannon$data %>% rownames_to_column("sample.id") # this moves the sample names to a new column that matches the metadata and allows them to be merged
#shows that all 39 samples are used. 
gplots::venn(list(metadata=metadata$sample.id, shannon=shannon$sample.id))

#need to separate labels for shannon
new_list = list()
samples=list(shannon$sample.id)
  
for (val in samples){
  new_val = str_sub(val, -0, - 2)
  new_list <- new_val
}
shannon$sample.type <- new_list

# Get mean & standard deviation by group
#displays a dataframe with sampletyp mean and sd
data_shannon <- shannon %>%
  group_by(sample.type) %>%
  summarise_at(vars(shannon_entropy),
               list(mean = mean,
                    sd = sd)) %>% 
  as.data.frame()
data_shannon  

# bargraph with error bars 
#ggplot(data_shannon) +
#  geom_bar( aes(x=sample.type, y=mean), stat="identity") +
#  geom_errorbar( aes(x=sample.type, y=mean, ymin=mean-sd, ymax=mean+sd), colour="black")

#gg boxplot with outlier
# Change outlier, color, shape and size
ggplot(shannon, aes(x=sample.type, y=shannon_entropy)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

# Compute the analysis of variance ANOVA
#We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest. 
res.aov <- aov(shannon_entropy ~ sample.type, data = shannon)
# Summary of the analysis
summary(res.aov)
#That tells us that there is a difference, but does not tell us which means are different. 
#It is a multiple range test similar to the LSD test except that Tukey utilized the honestly significant difference (HSD) test or the w-procedure.
#Any two treatments that mean having a difference more than honestly significant difference are said to be significantly different, otherwise not. 
TukeyHSD(res.aov, conf.level = .95)
#plot of HSD for groups - visualization of confidence intervals. 
par(mar = c(20, 8, 4, 2) + 0.1)  
plot(TukeyHSD(res.aov, conf.level=.95), las = 2)

#this prints out summary statisitics
shannon %>% 
  group_by(sample.type) %>%
  get_summary_stats(shannon_entropy, type = "common")
#this displays a box plot 
ggboxplot(shannon, x = "sample.type", y = "shannon_entropy")
#the displays non-parametric KW Test across all samples
res.kruskal <- shannon %>% kruskal_test(shannon_entropy ~ sample.type)
res.kruskal
#ensures large enough effect size to use KW 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
shannon %>% kruskal_effsize(shannon_entropy ~ sample.type)
#pairwise comparisons from KW with Dunn
# dunn test using bonferroni correction method
pwc <- shannon %>% 
  dunn_test(shannon_entropy ~ sample.type, p.adjust.method = "bonferroni") 
pwc
#pairwise comparison from KW with wilcox
pwc2 <- shannon %>% 
  wilcox_test(shannon_entropy ~ sample.type, p.adjust.method = "bonferroni")
pwc2
#visualize box plot with pairwise comparison
# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "sample.type")
ggboxplot(shannon, x = "sample.type", y = "shannon_entropy") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#evenness statistical test
#evenness is an alpha diversity metric that measures distribution of abundances of the group
evenness<-read_qza("/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/8_23_2022_ROutputs/core_metrics/evenness_vector.qza")
evenness_data = evenness$data
#need to separate labels for evenness
new_list = list()
samples=list(rownames(evenness_data))

for (val in samples){
  new_val = str_sub(val, -0, - 2)
  new_list <- new_val
}
evenness_data$sample.type <- new_list

#gg boxplot with outlier
# Change outlier, color, shape and size
ggplot(evenness_data, aes(x=sample.type, y=pielou_evenness)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

#this prints out summary statisitics
evenness_data %>% 
  group_by(sample.type) %>%
  get_summary_stats(pielou_evenness, type = "common")
#this displays a box plot 
ggboxplot(evenness_data, x = "sample.type", y = "pielou_evenness")
#the displays non-parametric KW Test across all samples
res.kruskal <- evenness_data %>% kruskal_test(pielou_evenness ~ sample.type)
res.kruskal
#ensures large enough effect size to use KW 
#The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
evenness_data %>% kruskal_effsize(pielou_evenness ~ sample.type)
#pairwise comparisons from KW with Dunn
# dunn test using bonferroni correction method
pwc <- evenness_data %>% 
  dunn_test(pielou_evenness ~ sample.type, p.adjust.method = "bonferroni") 
pwc
#pairwise comparison from KW with wilcox
pwc2 <- evenness_data %>% 
  wilcox_test(pielou_evenness ~ sample.type, p.adjust.method = "bonferroni")
pwc2
#visualize box plot with pairwise comparison
# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "sample.type")
ggboxplot(evenness_data, x = "sample.type", y = "pielou_evenness") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


##functions fun
boxplot_outlier <- function(df, df_dataCol, df_labelCol) {
  #gg boxplot with outlier
  # Change outlier, color, shape and size
  box_outlier = ggplot(df, aes(x=df_labelCol, y=df_dataCol)) + 
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4)
  ggsave(path = "/Users/suziepalmer/Desktop/Gut_BSI_16S_Data/alpha/"+, device='tiff', dpi=700)
  return(box_outlier)
}

boxplot_outlier(evenness_data, evenness_data$pielou_evenness, evenness_data$sample.type, )

