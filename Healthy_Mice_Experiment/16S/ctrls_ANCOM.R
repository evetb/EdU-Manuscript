# Purpose ---------------------------------------------------------------------------------------------

# ANCOM analysis for the "healthy mice" experiment in the EdU manuscript

# Libraries ------------------------------------------------------------------------------------------- 

library(phyloseq)
library(nlme)
library(compositions)
library(tidyverse)
library(qiime2R)

# QIIME2R ---------------------------------------------------------------------------------------------

# Loading the QIIME output files into phyloseq using QIIME2R

#### Healthy mice (controls) ####

ctrls_physeq <- qza_to_phyloseq(features = "ns-table-no-contam.qza",
                                    tree = "rooted-tree.qza",
                                    taxonomy = "taxonomy.qza",
                                    metadata = "ctrls_metadata.tsv")

ctrls_meta <- read.table("ctrls_metadata.tsv", header = TRUE)

# ANCOM -----------------------------------------------------------------------------

# Conducting differential abundance analysis (DAA) on the results from the healthy mice experiment
# We want to see whether there are differentially abundant species between
# the different sorted fractions (EdU+, EdU-, Whole)

#### Preparing ANCOM input ####

# Function to get filtered feature table based on specific taxonomic level (to input into ANCOM)

# Takes the phyloseq object that you want to use for filtering
# (ie - whole experiment, or just a specific time period/treatment)
# and a string indicating the taxonomic level of interest

ft_tax_filter <- function(phylo, tax_level) { 
  
  tax_filter <- tax_glom(phylo, tax_level, NArm = FALSE) 
  tax_filter_ft <- as(otu_table(tax_filter), "matrix")
  
}

# Filtering the phyloseq object to different
# Taxonomic levels

phylum_ctrls <- ft_tax_filter(ctrls_physeq, "Phylum")
fam_ctrls <- ft_tax_filter(ctrls_physeq, "Family")
genus_ctrls <-  ft_tax_filter(ctrls_physeq, "Genus")
species_ctrls <-  ft_tax_filter(ctrls_physeq, "Species")

#### ANCOM Function #### 

# Function to use ANCOM within each sorted fraction (sf), 
# accounting for experiment (abx1, abx2) as a random effect

ANCOM_sf_exp <- function(meta, ft, experiment_ft, csv_name) {
  source("ancom_v2.1.R")
  meta_data = meta
  feature_table = ft
  sample_var = "sample.id"
  group_var = NULL
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, 
                                     group_var, out_cut, zero_cut, lib_cut, neg_lb)
  
  feature_table = prepro$feature_table # Pre-processed feature table
  meta_data = prepro$meta_data # Pre-processed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  lme_control = NULL
  main_var = 'sf'
  p_adj_method = "BH"
  alpha = 0.05
  adj_formula = NULL
  rand_formula = '~ 1| experiment' 
  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
              alpha, adj_formula, rand_formula, lme_control)
  
  res$out <- arrange(res$out, desc(detected_0.6))
  
  taxa_list <- as.data.frame(tax_table(experiment_ft))
  
  taxa_list <- rownames_to_column(taxa_list, var = "taxa_id") 
  
  taxa_output <- left_join(res$out, taxa_list, by = "taxa_id")
  
  write.csv(taxa_output, csv_name)
  
}

ANCOM_sf_exp(ctrls_meta, phylum_ctrls, ctrls_physeq, "ANCOM_ctrls_sf_phy.csv")
ANCOM_sf_exp(ctrls_meta, fam_ctrls, ctrls_physeq, "ANCOM_ctrls_sf_fam.csv")
ANCOM_sf_exp(ctrls_meta, genus_ctrls, ctrls_physeq, "ANCOM_ctrls_sf_gen.csv")
ANCOM_sf_exp(ctrls_meta, species_ctrls, ctrls_physeq, "ANCOM_ctrls_sf_sp.csv")