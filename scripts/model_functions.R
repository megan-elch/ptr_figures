# function to compute necessary quantities to run model
compute_suff_stats = function(df, df_info, type = "mRNA", clusters = NULL, pop_label = NULL, already_logged = T, protein_label = "prot", peptide_label = "pep", omit_zeros = FALSE){
  df_info = na.omit(df_info)
  # select modality
  if(tolower(type) == "mrna"){
    df = df %>%
      as.matrix() %>%
      data.frame()
    if(omit_zeros == TRUE){
      df = replace(df, df==0, NA) # 0s treated as NAs
    }
    df$SYMBOL = rownames(df)
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id") %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      base::merge(df_info) %>% # merge cell type labels
      filter(ct %in% clusters) %>% # filter observations to be in correct clusters
      dplyr::group_by(id) %>%
      dplyr::mutate(ts_counts = sum(value[is.finite(value)])) %>% # total number of transcripts in each cell
      ungroup() %>%
      group_by(SYMBOL, ct) %>%
      dplyr::summarise(n_cells = sum(is.finite(value)),
                       mrna_sum = sum(value, na.rm = T), # models sum of transcripts counts for each gene, cell type
                       mrna_av = mean(value, na.rm = T),
                       mrna_mean_log2 = mean(log2(value + 1), na.rm = T),
                       counts = sum(ts_counts)) %>% # total number of transcript counts across cells in each cluster
      ungroup()

    if(!is.null(pop_label)){
      means = means %>% mutate(pop_mrna = pop_label)
    }
  }

  if(tolower(type) == "protein"){
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id")
    if(already_logged == FALSE){ # take log2 if needed
      means = means %>% mutate(value = log2(value))
    }
    means[["UNIPROT"]] = means[[protein_label]]
    means[["pep"]] = means[[peptide_label]]
    means = means %>%
      base::merge(df_info) %>%
      dplyr::group_by(UNIPROT, ct) %>%
      dplyr::mutate(n_cells_prot = length(unique(id[is.finite(value)]))) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      dplyr::group_by(pep, ct) %>%
      dplyr::summarise(pep_av = mean(value, na.rm = T), # models average log2 peptide intensity (original data already log2 transformed)
                       n_cells = sum(is.finite(value)), # number of cells for each peptide, cell type
                       n_cells_prot = unique(n_cells_prot), # number of cells for each protein, cell type
                       UNIPROT = unique(as.character(UNIPROT)),
                       pep_sum = sum(value, na.rm = T),
                       pep_sum_sq = sum(value^2, na.rm = T)) %>% # sum and sum sq if we want to increment target likelihood directly
      filter(n_cells > 5, ct %in% clusters) %>% # filter observations to be in correct clusters
      dplyr::group_by(UNIPROT) %>%
      dplyr::mutate(npep_ct = length(unique(pep))) %>% # number of peptides observed across cell types
      ungroup()

    if(!is.null(pop_label)){
      means = means %>% mutate(pop_protein = pop_label)
    }
  }
  return(means)
}

# prepare data to be inputted into model
model_prep = function(mrna_suff, protein_suff,
                      n_mrna = 2, n_protein = 2, min_ct = 3,
                      clusters, specific_genes = NULL){
  options(mc.cores = parallel::detectCores())

  protein_suff = protein_suff %>% mutate(pop_factor_protein = as.numeric(factor(pop_protein, levels = paste0("Pop", 1:n_protein))))
  mrna_suff = mrna_suff %>% mutate(pop_factor_mrna = as.numeric(factor(pop_mrna, levels = paste0("Pop", 1:n_mrna))))

  # filter out genes with no observed transcripts for all cell types
  mrna_suff = na.omit(mrna_suff) %>%
    dplyr::group_by(SYMBOL, pop_mrna)  %>%
    dplyr::mutate(n_ct = length(unique(ct)), zero_all = sum(mrna_sum) == 0) %>%
    ungroup() %>%
    filter(n_ct >= min_ct & zero_all == FALSE) %>%
    dplyr::select(-c(n_ct, zero_all))

  mrna_cts = mrna_suff %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(ct_m = unique(ct)) %>% # keeping track of number of cell types
    ungroup()

  protein_suff = na.omit(protein_suff) %>%
    dplyr::group_by(UNIPROT, pop_protein) %>%
    dplyr::mutate(n_ct = length(unique(ct)),
                  var_av = var(pep_av, na.rm = T)) %>% # compute number of cell types observed for each gene, data set
    ungroup() %>%
    dplyr::filter(n_ct >= min_ct & is.finite(log(var_av))) %>% # must have enough cell types observed
    dplyr::select(-n_ct)

  print(head(protein_suff))

  protein_cts = protein_suff %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(ct_p = unique(ct)) %>% # keeping track of number of cell types
    ungroup()

  cts = merge(mrna_cts, protein_cts) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(n_ct_int = length(intersect(unique(ct_p), unique(ct_m)))) %>% # compute number of cell types shared across modalities
    ungroup() %>%
    dplyr::filter(n_ct_int >= min_ct) %>% # keep only genes observed in enough cell types
    pull(UNIPROT)

  celltypes_union = intersect(mrna_suff$ct, protein_suff$ct) # keep shared cell types

  proteins_union = protein_suff %>%
    filter(UNIPROT %in% mrna_suff$UNIPROT, UNIPROT %in% cts) %>% # keep shared proteins in enough cell types
    pull(UNIPROT)

  print(length(proteins_union))

  if(!is.null(specific_genes)){
    proteins_union = intersect(proteins_union, specific_genes)
  }

  UNIPROT_map = mrna_suff %>% # make reference df to keep track of gene identities in model
    dplyr::filter(UNIPROT %in% proteins_union, ct %in% celltypes_union, UNIPROT %in% cts) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise() %>%
    ungroup() %>%
    dplyr::mutate(uni_map = 1:length(unique(UNIPROT))) # id number to go with each gene

  print(paste("number of genes in model:", nrow(UNIPROT_map)))

  mrna_suff = mrna_suff %>% filter(ct %in% celltypes_union) # cell type filtering
  protein_suff = protein_suff %>% filter(ct %in% celltypes_union) # cell type filtering

  # celltype mapping, population mapping
  protein_suff_sample = merge(protein_suff, UNIPROT_map) %>% # merge to remove all filtered out obs
    arrange(pop_protein) %>%
    mutate(ct_factor = as.numeric(factor(ct, levels = clusters))) # make factor for cell type ids
  print(paste("Number of Peptide Observations:", nrow(protein_suff_sample)))

  mrna_suff_sample = merge(mrna_suff, UNIPROT_map) %>%
    arrange(pop_mrna) %>%
    mutate(ct_factor = as.numeric(factor(ct, levels = clusters)))

  print(paste("Number of mRNA Observations:", nrow(mrna_suff_sample)))
  list(mrna = mrna_suff_sample, protein = protein_suff_sample, gene_map = UNIPROT_map, clusters = clusters)
}

# run stan model
run_model = function(prep_list, project_name, n_iter = 2000, n_warmup = 1000, n_chains = 5,
                     seed_id = 100, n_cores = 5, refresh = 100, init = 0){

  mrna = prep_list$mrna
  protein = prep_list$protein
  gene_map = prep_list$gene_map
  clusters = prep_list$clusters
  total_obs_mrna = nrow(mrna)
  total_obs_protein = nrow(protein)

  stan_input = list(n_g = nrow(gene_map), # number of genes
                    total_obs_mrna = total_obs_mrna, # total number of observations (mrna)
                    total_obs_protein = total_obs_protein, # total number of observations (protein)
                    mrna_sum = mrna$mrna_sum, # observed mrna transcript sum
                    protein_av = protein$pep_av, # observed log2 peptide averages
                    mrna_gene_rec = as.numeric(mrna$uni_map), # keep track of gene associated w each observation
                    protein_gene_rec = as.numeric(protein$uni_map),  # keep track of protein associated w each observation
                    mrna_celltype_rec = as.numeric(mrna$ct_factor), # keep track of cell types mrna
                    protein_celltype_rec = as.numeric(protein$ct_factor), # keep track of cell types protein
                    mrna_pop_rec = as.numeric(mrna$pop_factor_mrna), # keep track of mrna data set
                    protein_pop_rec = as.numeric(protein$pop_factor_protein), # keep track of protein data sets
                    n_mrna = mrna$counts, # total transcript counts
                    n_protein = protein$n_cells) # number of cells protein

  print("sampling starting...")
  ptm = proc.time() # run model
  stan_fit = rstan::stan(file=paste0("scripts/", project_name, "/", project_name, ".stan"),
                         data=stan_input, iter = n_iter, warmup = n_warmup,
                         chains=n_chains, seed=seed_id, cores = n_cores,
                         refresh = refresh, init = init)
  print(proc.time() - ptm)

  list(stan_fit = stan_fit, gene_map = gene_map, clusters = clusters)
}
