library(tidyverse)
library(Matrix)
library(rlist)
library(readxl)
library(org.Hs.eg.db)  

# identify project
project_name = "testes"
source("scripts/model_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St") # identify cell types order

# load previously filtered mrna data
load(paste0("data/scRNA_seq/mrna_noSeurat.RData"))

# load in peptide data sets 
peptides_pop1 = read.table(paste0("data/scProtein/Pop1/testis_nPoP.181122_peptideLevelMatrix_dartUpdate_medIntNorm_BCAtt3.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop2 = read.table(paste0("data/scProtein/Pop2/testis_nPoP.081221_peptideLevelMatrix_dartUpdate.medIntNorm.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop3 = read.table(paste0("data/scProtein/Pop3/testis_nPoP.130122_peptideLevelMatrix_dartUpdate_medIntNormNoBC.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop4 = read.table(paste0("data/scProtein/Pop4/testis_nPoP.030222_peptideLevelMatrix_dartUpdate_medIntNorm.txt"), 
                           sep = '\t', header = TRUE)

# protein meta data
protein_pop1_meta = read.delim(paste0("data/scProtein/Pop1/liger.clustToCellID_N4Sep.1_k10.L5_dartUpdateNewNormNoBC_t.1.4.noProteo.txt"))
rownames(protein_pop1_meta) = protein_pop1_meta$id
protein_pop2_meta = read.delim(paste0("data/scProtein/Pop2/liger.clustToCellID_N4Sep.2_k10.L5_dartUpdateNewNormNoBC_t1.4.noProteo.txt"))
rownames(protein_pop2_meta) = protein_pop2_meta$id
protein_pop3_meta = read.delim(paste0("data/scProtein/Pop3/liger.clustToCellID_N4Sep.3_k10.L5_dartUpdateNewNormNoBC_protAllRows.noProteo.txt"))
rownames(protein_pop3_meta) = protein_pop3_meta$id
protein_pop4_meta = read.delim(paste0("data/scProtein/Pop4/liger.clustToCellID_N4Sep.4_k10.L45_newNorm_t1.4.PlusClustMarks.txt"))
rownames(protein_pop4_meta) = protein_pop4_meta$id

# reformat protein meta data
protein_pop1_info = data.frame(celltype = protein_pop1_meta$cellType, 
                               id = rownames(protein_pop1_meta))
protein_pop2_info = data.frame(celltype = protein_pop2_meta$cellType, 
                               id = rownames(protein_pop2_meta))
protein_pop3_info = data.frame(celltype = protein_pop3_meta$cellType, 
                               id = rownames(protein_pop3_meta))
protein_pop4_info = data.frame(celltype = protein_pop4_meta$cellType, 
                               id = rownames(protein_pop4_meta))
protein_pop1_info$ct = protein_pop1_info$celltype
protein_pop2_info$ct = protein_pop2_info$celltype
protein_pop3_info$ct = protein_pop3_info$celltype
protein_pop4_info$ct = protein_pop4_info$celltype

p = mrna_pop1_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

p = mrna_pop2_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

p = protein_pop1_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

p = protein_pop2_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

p = protein_pop3_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

p = protein_pop4_info %>% dplyr::group_by(ct) %>% dplyr::summarise(ct_sum = length(unique(id)))
print(p)

# compute sufficient statistics for each data set
mrna_suff_pop1 = compute_suff_stats(mrna_pop1_counts, mrna_pop1_info, clusters = clusters, pop_label = "Pop1") 
print("mrna1")
mrna_suff_pop2 = compute_suff_stats(mrna_pop2_counts, mrna_pop2_info, clusters = clusters, pop_label = "Pop2") 
print("mrna2")
mrna_suff = rbind(mrna_suff_pop1, mrna_suff_pop2)
print(head(mrna_suff))

protein_suff_pop1 = compute_suff_stats(peptides_pop1, protein_pop1_info, type = "protein", clusters = clusters, pop_label = "Pop1") 
print("protein1")
protein_suff_pop2 = compute_suff_stats(peptides_pop2, protein_pop2_info, type = "protein", clusters = clusters, pop_label = "Pop2") 
print("protein2")
protein_suff_pop3 = compute_suff_stats(peptides_pop3, protein_pop3_info, type = "protein", clusters = clusters, pop_label = "Pop3") 
print("protein3")
protein_suff_pop4 = compute_suff_stats(peptides_pop4, protein_pop4_info, type = "protein", clusters = clusters, pop_label = "Pop4") 
print("protein4")
protein_suff = rbind(protein_suff_pop1, protein_suff_pop2, 
                     protein_suff_pop3, protein_suff_pop4)
print(head(protein_suff))
 
# link together labels for each modality
protein_genes_union = unique(protein_suff$UNIPROT)
mrna_genes_union = unique(mrna_suff$SYMBOL)

protein_genes = AnnotationDbi::select(org.Hs.eg.db, keys = protein_genes_union, 
                                      keytype ='UNIPROT', columns = c("SYMBOL"))
protein_genes = protein_genes[!duplicated(protein_genes[,"UNIPROT"]),]
mrna_genes = AnnotationDbi::select(org.Hs.eg.db, keys= mrna_genes_union, 
                                   keytype = 'SYMBOL', columns = c("UNIPROT"))
 
# combine with sufficient statistics                                       
all_genes = merge(protein_genes, mrna_genes)
protein_suff = merge(protein_suff, all_genes)
mrna_suff = merge(mrna_suff, all_genes)

save(mrna_suff, protein_suff, file = paste0("processed_data/suff_stats_", project_name, ".RData"))