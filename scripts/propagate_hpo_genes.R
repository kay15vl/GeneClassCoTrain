#Created by: Karla Vela Lopez
#Created on: 4/1/24

# Input: HPO data
# Output: Propagated genes associated to phenotype of diseases of interest 
# (TB phenotypes).

library(tidyverse)
library(ontologyIndex)
library(parallel)
library(data.table)

#read in ontology  files
hpo_ont <- get_OBO(file = "data/raw/hp.obo",
                   propagate_relationships = "is_a",
                   extract_tags = "everything")
df_hpo_ont <- as.data.frame(hpo_ont)

mondo_ont <- get_OBO(file = "data/raw/mondo_2024_03.obo",
                     propagate_relationships = "is_a",
                     extract_tags = "everything")
df_mondo_ont <- as.data.frame(mondo_ont)

#hpo annotations file
hpoa_file_path <- "data/raw/phenotype.hpoa"
hpoa_df <- read.delim(hpoa_file_path, header = TRUE, 
                      skip = 4, stringsAsFactors = FALSE) %>% 
  rename(disease_id = database_id) %>% 
  select(-reference:-biocuration, -qualifier) %>% setDT()

# phenotypes not including ancestors
hpgenestopheno <-"data/raw/genes_to_phenotype.txt"
hp_genes_to_pheno <-read_delim(hpgenestopheno,
                                                      delim = "\t",
                               escape_double = FALSE,
                                                      trim_ws = TRUE) %>% 
  select(-frequency)  %>% setDT()



map_mondo_to_ontology <- function(dis_gene_data, obo_file) {

  mondo_to_disid <- as.data.frame(obo_file) %>%
    dplyr::select(
      -`format-version`:-ontology, -synonym, -comment:-union_of,
      -seeAlso:-is_class_level, -parents:-ancestors, -property_value:-def,
      -replaced_by:-is_a
    ) %>%
    mutate(xref = lapply(xref, function(x) gsub("; ", ";", x))) %>%
    separate_rows(xref, sep = ";") %>%
    filter(!is.na(obsolete) & obsolete != TRUE) %>%
    dplyr::rename(mondo_id = id, mondo_term = name, ontology_dis_id = xref) %>%
    filter(!is.na(ontology_dis_id) & ontology_dis_id != "") %>% 
    rename(disease_id = ontology_dis_id) %>% 
    mutate(disease_id = str_replace(disease_id, "Orphanet", "ORPHA"))
  
  # join mondo df to the gwas dataframe. Join used to map mondo to the dise
  convert_to_mondo <- dis_gene_data %>%
    left_join(mondo_to_disid, by = "disease_id") %>%
    relocate(c("mondo_id", "mondo_term"), .after = disease_id) %>%
    dplyr::select(-obsolete, -subset) %>%
    filter(!is.na(mondo_id)) %>%
    unique() 
  
  return(convert_to_mondo)
}

#funciton to get ancestor information
get_ances<- function(term_id, ontology_obo){
  
  # get parent terms
  ancestors<- get_ancestors(ontology_obo, terms= term_id )
  
  # init df
  parent_df <- tibble()
  
  # add data to df
  parent_df <- tibble(ID= term_id,
                      Parent= ancestors)
  
  
  # collapse ancestors to one row per go term id
  parent_df <- parent_df %>% filter(Parent != term_id) %>% 
    group_by(ID) %>% 
    summarise(parent_ct = n_distinct(Parent),
              Parent = paste(Parent, collapse = ","))
  return(parent_df)
  
}

# Function to propagate genes from child terms to parent terms
propagate_genes<- function(current_term, ance_df, disgenes_df, df_mondo_ont){
  
  # get the genes of the current term
  current_term_genes <- disgenes_df %>% filter(mondo_id == current_term)
  
  #get ancestor terms
  ancestor_terms <- ance_df %>% filter(ID == current_term) %>%select(Parent) %>%
    unlist() %>% strsplit( ",") %>% unlist() %>%trimws() %>%  as.vector()
  
  # check to see which ancestors are in the hpo disease df
  ancestors_in_hpo_disdf <- disgenes_df %>% filter(mondo_id %in% ancestor_terms)
  
  
  # check if there are ancestors if not make propagated genes from scratch
  if(nrow(ancestors_in_hpo_disdf)==0){
    get_ancestors_from_ontology <- df_mondo_ont %>% filter(id %in% ancestor_terms)
    ancestors <- data.frame(
      mondo_id = rep(unique(get_ancestors_from_ontology$id), each = nrow(current_term_genes)),
      mondo_term = rep(unique(get_ancestors_from_ontology$name), each = nrow(current_term_genes)),
      gene_id = rep(current_term_genes$gene_id),
      hpo_id = rep(current_term_genes$hpo_id)
    )
    
    # add to the rest of the ancestor data
    propagated_genes <- bind_rows(current_term_genes, ancestors)
    
  }
  
  # get the missing ancestor terms if any
  missing_ancestors <- setdiff(ancestor_terms, ancestors_in_hpo_disdf$mondo_id)
  
  if(length(missing_ancestors)>0 & nrow(ancestors_in_hpo_disdf)!=0){
    #create data frame of missing ancestors
    
    # get the id and name of the missing mondo ancestors
    get_missing_ancestor_identifiers <- df_mondo_ont %>% filter(id %in% missing_ancestors)
    
    missing_ancestors <- data.frame(
      mondo_id = rep(unique(get_missing_ancestor_identifiers$id), each = nrow(current_term_genes)),
      mondo_term = rep(unique(get_missing_ancestor_identifiers$name), each = nrow(current_term_genes)),
      gene_id = rep(current_term_genes$gene_id),
      hpo_id = rep(current_term_genes$hpo_id)
    )
    
    # propagated genes to ancestors
    propagate <- data.frame(
      mondo_id = rep(unique(ancestors_in_hpo_disdf$mondo_id), each = nrow(current_term_genes)),
      mondo_term = rep(unique(ancestors_in_hpo_disdf$mondo_term), each = nrow(current_term_genes)),
      gene_id = rep(current_term_genes$gene_id),
      hpo_id = rep(current_term_genes$hpo_id)
    )
    
    # add to the rest of the ancestor data
    propagated_genes <- bind_rows(missing_ancestors, propagate, 
                                  ancestors_in_hpo_disdf, 
                                  current_term_genes)
  }
  
  
  # keep unique rows
  final_propagated_genes <- propagated_genes %>% 
    distinct(mondo_id, gene_id, hpo_id, .keep_all = TRUE)
  
  return(final_propagated_genes)
}



####  HPO Gene Collection ####

# get phenotype data
complete_hpo <- hpoa_df %>% inner_join(hp_genes_to_pheno,
                                      by=c("disease_id", "hpo_id")) %>% 
  dplyr::rename("gene_id"= "ncbi_gene_id") %>% 
  select(-gene_symbol) %>% relocate(gene_id, .after ="disease_name") %>% 
  unique()

# separate the disease data and phenotype data
disease_info_hpo <- complete_hpo %>% select(disease_id:hpo_id) %>% unique()

# get ancestor information for the diseases
# map mondo identifiers to the disease ids in hpo
get_mondo_hpo <- map_mondo_to_ontology(disease_info_hpo, df_mondo_ont)  %>% 
  distinct(mondo_id, gene_id, hpo_id, .keep_all = TRUE) %>% select(-disease_name, -disease_id)

# get the terms to propagate them
hpo_mondo_terms <- unique(get_mondo_hpo$mondo_id)

#get ancestor information
hpo_mondo_ances <- mclapply(hpo_mondo_terms, get_ances, mondo_ont) %>% bind_rows() %>% 
  unique() %>% setDT()

# propagate the disease terms
propagate_mondo_hpo <-  mclapply(hpo_mondo_terms, propagate_genes, hpo_mondo_ances, get_mondo_hpo, df_mondo_ont) %>% 
  bind_rows() %>%  distinct(mondo_id, gene_id, hpo_id, .keep_all = TRUE) %>% 
  setDT()

rm(hpoa_df, disease_info_hpo, hpo_mondo_ances)
gc()


# get tb phenotypes HP:0032262- pulmonary tb, HP:6000542- Positive CSF mycobacterium tuberculosis nucleic acid test
# MONDO:0019146- inherited susceptibility to mycobacterial diseases
# HP:0032271 - extrapulmonary tuberculosis
tb_hpos <- propagate_mondo_hpo %>% 
  filter(mondo_id %in% c("MONDO:0030491","MONDO:0019146")) %>% 
  inner_join(complete_hpo, by= "hpo_id") %>% 
  select(-disease_id:-gene_id.y) %>% 
  rename(gene_id = gene_id.x) %>% 
  distinct(gene_id, hpo_id, .keep_all = TRUE)

# since the gene list is so short, I decided to get more genes that are involved with
# the hpos that were in the tb_hpos set. Since, these hpos fall under the mondo
# ids of interest, I gathered the other genes associated to the hpos. 
tb_other_hpos <- unique(tb_hpos$hpo_id)
more_hpo_genes <- propagate_mondo_hpo %>% filter(hpo_id %in% tb_other_hpos) %>% 
  inner_join(complete_hpo, by = "hpo_id") %>% 
  select(-disease_id:-gene_id.y) %>% 
  rename(gene_id = gene_id.x) %>% 
  distinct(gene_id, hpo_id, .keep_all = TRUE)

#combine other hpo genes to the tb_hpos
final_tb_hpos <- bind_rows(tb_hpos, more_hpo_genes)

save(final_tb_hpos, file = "data/processed/propagated_hpo_disgenes.RData")

# obtain only the genes and gene ids
propagated_hpo_genes <- final_tb_hpos %>% 
  distinct(gene_id, mondo_id, .keep_all = TRUE) %>% 
  select(hpo_id, hpo_name, gene_id) 
  
write_tsv(propagated_hpo_genes, file= "data/processed/propagated_hpo_genes.tsv")


print("Script successfully executed")








