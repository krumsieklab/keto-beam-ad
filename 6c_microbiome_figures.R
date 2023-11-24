# Script to generate figure panels  2c from the paper. 
# Generate pdfs with plots in the results folder

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries
library(readxl)
library(maplet)
library(dplyr)
library(MetBrewer)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(Hmisc)
library(Rmisc)
library(igraph)
library(visNetwork)
library(stringr)
source("custom_functions.R") # customized functions
# get colors
palette <- c(MetBrewer::met.brewer("Archambault"), "#0c7156")


# microbiome results
association_df <- read.csv("results/tmp_plot_one_bcaa_gene_encoded.csv")

##### Figure 2c: Microbiome plot #####


# Get color and label of microbiome species
association_df %<>% 
  dplyr::mutate(label=ifelse(subset_color_spp, species,NA)) %>% 
  dplyr::mutate(size = ifelse(subset_color_spp, 7,5))

colors <- c("grey", palette)
names(colors)<- association_df$label %>% 
  unique()

# Figure 2ci microbiome-BCAA scatterplot
fig_2ci_plot <- ggplot(association_df, aes(x= bcaa, y = diet_nocross.T.Ketogenic.Diet._mean, color=label, size = size)) + 
  
  geom_point()+
  
  theme_minimal()+
  
  scale_color_manual(values=colors)+
  
  scale_size_continuous(range = c(2,4))+
  
  xlab("Association to plasma BCAA")+
  
  ylab("Log-fold enrichment in MMKD over AHAD") +
  
  theme(legend.text = element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))+
  
  guides(colour = guide_legend(override.aes = list(size=8)))+
  
  ylim(-3,1.5)

fig_2c_plot_df <- melt(association_df) %>% 
  dplyr::filter(variable != "size") 

# Figure 2cii Box plot of species that encode BCAA biosynthesis
fig_2cii_plot <- ggplot(fig_2c_plot_df, aes(x=all_genes_bool, y = value, fill = all_genes_bool))+

  geom_boxplot()+
  
  scale_fill_manual(name = "Genome Encodes BCAA biosynthesis", 
                    labels = c("No","Yes"),
                    values=c(palette[2], palette[4])) +
  
  geom_signif(
    comparisons = list(c("No", "Yes")),
    annotations = "*",
    textsize = 10,
    vjust=0.5
  ) +
  
  facet_wrap(~variable,scales="free")+
  
  theme_bw() +
  
  xlab("") +
  
  ylab("Association to Plasma BCAA")+
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        strip.text = element_blank(),
        strip.background = element_rect(fill="seashell"))

 pdf("results/Figure_2ci.pdf", height=4, width=8)
 print(fig_2ci_plot)
 dev.off()

pdf("results/Figure_2cii.pdf")
print(fig_2cii_plot)
dev.off()
