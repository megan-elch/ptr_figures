# function to compute necessary quantities to run model
compute_suff_stats = function(df, df_info, type = "mRNA", clusters = NULL, pop_label = NULL, already_logged = T, protein_label = "prot", peptide_label = "pep", omit_zeros = FALSE){
  df_info = na.omit(df_info)
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
      dplyr::summarise(n_cells = sum(value > 0),
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

compute_reliability = function(df, df_info, protein_label = "prot", peptide_label = "pep", seed_id = 100){
    df[["UNIPROT"]] = df[[protein_label]]
    df[["pep"]] = df[[peptide_label]]
    cells = intersect(df_info$id, colnames(df))
    
    set.seed(seed_id)
    df_group_a = df %>% 
        dplyr::group_by(UNIPROT) %>%
        dplyr::sample_n(size = floor(length(unique(pep))/2)) %>%
        ungroup()
        
    peptides_a = unique(df_group_a$pep)
    peptides_b = setdiff(unique(df$pep), peptides_a)
    
    df_group_a = df_group_a %>%
        pivot_longer(cols = cells, names_to = "id") %>%
        dplyr::group_by(UNIPROT, id) %>%
        dplyr::summarise(prot_a = mean(value, na.rm = T)) %>%
        ungroup()
        
    reliabilities = df %>%
        filter(pep %in% peptides_b) %>%
        pivot_longer(cols = cells, names_to = "id") %>%
        dplyr::group_by(UNIPROT, id) %>%
        dplyr::summarise(prot_b = mean(value, na.rm = T)) %>%
        ungroup() %>%
        merge(df_group_a) %>%
        dplyr::group_by(UNIPROT) %>%
        dplyr::summarise(reliability = cor(prot_a, prot_b, use = "pairwise.complete.obs"))
    return(reliabilities)
}

compute_reliability_ppc = function(protein_ppc, seed_id = 100){
    set.seed(seed_id)
    iter_a = unique(protein_ppc$.iteration) %>% sample(floor(length(unique(protein_ppc$.iteration)))/2)
    iter_b = setdiff(unique(protein_ppc$.iteration), iter_a)
    
    ppc_a = protein_ppc %>% 
        filter(.iteration %in% iter_a) %>%
        dplyr::group_by(UNIPROT, ct, pop_protein) %>%
        dplyr::summarise(prot_a = mean(protein_bar_rep, na.rm = T))
        
    ppc_b = protein_ppc %>%
        filter(.iteration %in% iter_b) %>%
        dplyr::group_by(UNIPROT, ct, pop_protein) %>%
        dplyr::summarise(prot_b = mean(protein_bar_rep, na.rm = T))
        
    ppc_reliabilities = ppc_a %>% merge(ppc_b) %>%
        dplyr::group_by(UNIPROT, pop_protein) %>%
        dplyr::summarise(reliability = cor(prot_a, prot_b, use = "pairwise.complete.obs"))
        
    return(ppc_reliabilities)
}

compute_reliability_ct = function(df, df_info, protein_label = "prot", peptide_label = "pep", seed_id = 100){
    df[["UNIPROT"]] = df[[protein_label]]
    df[["pep"]] = df[[peptide_label]]
    cells = intersect(df_info$id, colnames(df))
    
    set.seed(seed_id)
    df_group_a = df %>% 
        dplyr::group_by(UNIPROT) %>%
        dplyr::sample_n(size = floor(length(unique(pep))/2)) %>%
        ungroup()
        
    peptides_a = unique(df_group_a$pep)
    peptides_b = setdiff(unique(df$pep), peptides_a)
    
    df_group_a = df_group_a %>%
        pivot_longer(cols = cells, names_to = "id") %>%
        merge(df_info) %>%
        dplyr::group_by(UNIPROT, ct) %>%
        dplyr::summarise(prot_a = mean(value, na.rm = T)) %>%
        ungroup()
        
    reliabilities = df %>%
        filter(pep %in% peptides_b) %>%
        pivot_longer(cols = cells, names_to = "id") %>%
        merge(df_info) %>%
        dplyr::group_by(UNIPROT, ct) %>%
        dplyr::summarise(prot_b = mean(value, na.rm = T)) %>%
        ungroup() %>%
        merge(df_group_a) %>%
        dplyr::group_by(UNIPROT) %>%
        dplyr::summarise(reliability = cor(prot_a, prot_b, use = "pairwise.complete.obs"))
    return(reliabilities)
}

model_prep = function(mrna_suff, protein_suff, 
                      n_mrna = 2, n_protein = 2, min_ct = 3,
                      clusters, specific_genes = NULL){
    options(mc.cores = parallel::detectCores())
    
    protein_suff = protein_suff %>% 
        filter(ct %in% clusters) %>%
        mutate(pop_factor_protein = as.numeric(factor(pop_protein, levels = paste0("Pop", 1:n_protein))))
    mrna_suff = mrna_suff %>% 
        filter(ct %in% clusters) %>%
        mutate(pop_factor_mrna = as.numeric(factor(pop_mrna, levels = paste0("Pop", 1:n_mrna))))
    
    mrna_suff = na.omit(mrna_suff) %>%
                dplyr::group_by(SYMBOL, pop_mrna)  %>%
                dplyr::mutate(n_ct = length(unique(ct)), zero_all = sum(mrna_sum) == 0) %>%
                ungroup() %>%
                filter(n_ct >= min_ct & zero_all == FALSE) %>%
                dplyr::select(-c(n_ct, zero_all))
                
    print(head(mrna_suff))

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
                      
run_model = function(prep_list, project_name, n_iter = 2000, n_warmup = 1000, n_chains = 5,
                     seed_id = 100, n_cores = 5, refresh = 100, init = 0, n_cells = T){
    
    mrna = prep_list$mrna
    protein = prep_list$protein
    gene_map = prep_list$gene_map
    clusters = prep_list$clusters
    total_obs_mrna = nrow(mrna)
    total_obs_protein = nrow(protein)
    
    if(n_cells == FALSE){
        protein$n_cells = protein$npep_ct
    }
    
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

# function to return posterior samples of mu, r, and mu + r
extract_ratios = function(stan_fit, gene_map, clusters){
    posterior_draws = stan_fit %>%
                      spread_draws(mu[gene_num, celltype_num], r[gene_num, celltype_num]) %>%
                      ungroup() %>%
                      mutate(prot = mu + r,
                             UNIPROT = gene_map$UNIPROT[gene_num],
                             ct = clusters[celltype_num])
    return(posterior_draws)
}

# test r = 0 at gene level, using 95 pct interval or expected proportion of false discoveries
# put in posterior draws of mu, r, mu + r, or stan fit directly
# outputs data frame containing posterior summaries and gene, celltype test results
test_genes = function(posterior_draws = NULL, stan_fit = NULL, gene_map = NULL, clusters = NULL, threshold = 0.01){ 
    if(is.null(posterior_draws)){
        print("extracting ratios")
        posterior_draws = extract_ratios(stan_fit = stan_fit, gene_map = gene_map, clusters = clusters)
    }
    print(head(posterior_draws))
    gene_res = posterior_draws %>%
               na.omit() %>%
               dplyr::group_by(UNIPROT, ct) %>%
               dplyr::summarise(mu_av = mean(mu, na.rm = T), # posterior summaries mu
                                mu_lwr = quantile(mu, 0.025, na.rm = T),
                                mu_upr = quantile(mu, 0.975, na.rm = T),
                                mu_med = median(mu, na.rm = T),
                                r_av = mean(r, na.rm = T), # posterior summaries r
                                r_lwr = quantile(r, 0.025, na.rm = T),
                                r_upr = quantile(r, 0.975, na.rm = T),
                                r_med = median(r, na.rm = T),
                                prot_av = mean(prot, na.rm = T), # posterior samples prot
                                prot_lwr = quantile(prot, 0.025, na.rm = T),
                                prot_upr = quantile(prot, 0.975, na.rm = T),
                                prot_med = median(prot, na.rm = T),
                                significant = !(r_lwr <= 0 & r_upr >= 0)) %>% # interval size
                   ungroup() %>%
                   dplyr::group_by(UNIPROT) %>%
                   dplyr::mutate(n_ct = length(intersect(unique(ct[is.finite(mu_av)]), unique(ct[is.finite(prot_av)])))) %>%
                   ungroup()
    return(gene_res)
}

test_genes_centered = function(posterior_draws, threshold = 0.01){
    gene_res = posterior_draws %>%
               dplyr::group_by(UNIPROT, .iteration, .chain) %>%
               dplyr::mutate(r_gene = mean(r)) %>%
               ungroup() %>%
               dplyr::mutate(r = r - r_gene) %>%
               dplyr::group_by(UNIPROT, ct) %>%
               dplyr::summarise(p_r = ifelse(mean(r < 0) <= 1 - mean(r < 0), mean(r < 0), 1 - mean(r < 0)),
                                mu_av = mean(mu, na.rm = T),
                                mu_lwr = quantile(mu, 0.025, na.rm = T),
                                mu_upr = quantile(mu, 0.975, na.rm = T),
                                mu_med = median(mu, na.rm = T),
                                r_av = mean(r, na.rm = T), # posterior summaries r
                                r_lwr = quantile(r, 0.025, na.rm = T),
                                r_upr = quantile(r, 0.975, na.rm = T),
                                r_med = median(r, na.rm = T),
                                prot_av = mean(prot, na.rm = T), # posterior samples prot
                                prot_lwr = quantile(prot, 0.025, na.rm = T),
                                prot_upr = quantile(prot, 0.975, na.rm = T),
                                prot_med = median(prot, na.rm = T),
                                significant = !(r_lwr <= 0 & r_upr >= 0)) %>% # interval size
                   ungroup() %>%
                   dplyr::group_by(UNIPROT) %>%
                   dplyr::mutate(n_ct = length(intersect(unique(ct[is.finite(mu_av)]), unique(ct[is.finite(prot_av)])))) %>%
                   ungroup() %>%
                   dplyr::group_by(ct) %>%
                   arrange(p_r, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
                   dplyr::mutate(fdr = cummean(p_r[order(p_r)]), # expected proportion
                                 significant_fdr = fdr <= threshold) %>%
                   ungroup()
    return(gene_res)
}

test_genes_mean = function(posterior_draws, threshold = 0.01){
    gene_res = posterior_draws %>%
               na.omit() %>%
               dplyr::group_by(UNIPROT, .iteration, .chain) %>%
               dplyr::mutate(r_gene = mean(r)) %>%
               ungroup() %>%
               dplyr::group_by(UNIPROT, .iteration, .chain) %>% # record average across genes not associated with each group
               dplyr::mutate(r_anti = ((n()*r_gene) - r)/(n() - 1)) %>% 
               ungroup() %>%          
               dplyr::group_by(UNIPROT, ct) %>%
               dplyr::summarise(p_r = ifelse(mean(r < r_gene) <= 1 - mean(r < r_gene), mean(r < r_gene), 1 - mean(r < r_gene)),
                                p_r_anti = ifelse(mean(r < r_anti) <= 1 - mean(r < r_anti), mean(r < r_anti), 1 - mean(r < r_anti)),
                                mu_av = mean(mu, na.rm = T),
                                mu_lwr = quantile(mu, 0.025, na.rm = T),
                                mu_upr = quantile(mu, 0.975, na.rm = T),
                                mu_med = median(mu, na.rm = T),
                                r_av = mean(r, na.rm = T), # posterior summaries r
                                r_av_anti = mean(r_anti, na.rm = T),
                                r_av_gene = mean(r_gene, na.rm = T),
                                r_lwr = quantile(r, 0.025, na.rm = T),
                                r_upr = quantile(r, 0.975, na.rm = T),
                                r_med = median(r, na.rm = T),
                                prot_av = mean(prot, na.rm = T), # posterior samples prot
                                prot_lwr = quantile(prot, 0.025, na.rm = T),
                                prot_upr = quantile(prot, 0.975, na.rm = T),
                                prot_med = median(prot, na.rm = T)) %>% # interval size
                   ungroup() %>%
                   dplyr::group_by(UNIPROT) %>%
                   dplyr::mutate(n_ct = length(intersect(unique(ct[is.finite(mu_av)]), unique(ct[is.finite(prot_av)])))) %>%
                   ungroup() %>%
                   dplyr::group_by(ct) %>%
                   arrange(p_r, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
                   dplyr::mutate(fdr = cummean(p_r[order(p_r)]), # expected proportion
                                 significant = fdr <= threshold) %>%
                   ungroup() %>%
                   dplyr::group_by(ct) %>%
                   arrange(p_r_anti, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
                   dplyr::mutate(fdr_anti = cummean(p_r_anti[order(p_r_anti)]), # expected proportion
                                 significant_anti = fdr_anti <= threshold)                                 
    return(gene_res)
}

# return posterior means of scale parameter
# put in stan fit and gene mapping 
extract_scale = function(stan_fit, gene_map){
    scale_info = stan_fit %>% # use tidybayes to pull out scale parameter
                 spread_draws(scale_param[gene_num]) %>%
                 ungroup() %>%
                 mutate(UNIPROT = gene_map$UNIPROT[gene_num], #organize by gene
                        multiplier = 1/scale_param) # value that's actually multiplied by mu + r for each gene for get protein
    scale_means = scale_info %>% dplyr::group_by(UNIPROT) %>% dplyr::summarise(scale_av = mean(multiplier))
    
    return(scale_means)                    
}

# extract protein column normalization term. input stan fit and vector of cluster order
extract_b = function(stan_fit, clusters, n_protein = 2){
    protein_pop_order = paste0("Pop", 1:n_protein)
    b_info = stan_fit %>%
             spread_draws(b_protein[protein_pop_num, celltype_num]) %>%
             mutate(ct = clusters[celltype_num],
                    pop_protein = protein_pop_order[protein_pop_num])
    return(b_info)                
}

# function that returns posterior means of kappa for each gene, data set
extract_kappa = function(stan_fit, gene_map, n_protein){
    protein_pop_order = paste0("Pop", 1:n_protein)
    kappa_info = stan_fit %>%
                 spread_draws(kappa[protein_pop_num, gene_num]) %>% # extract posterior draws
                              mutate(UNIPROT = gene_map$UNIPROT[gene_num],
                                     pop_protein = protein_pop_order[protein_pop_num]) %>%
                 dplyr::group_by(pop_protein, UNIPROT) %>%
                 dplyr::summarise(kappa_av = mean(kappa)) # compute posterior mean
    return(kappa_info)
}

# function that returns posterior means of gamma for each gene, data set                   
extract_gamma = function(stan_fit, gene_map, n_mrna){
    mrna_pop_order = paste0("Pop", 1:n_mrna)
    gamma_info = stan_fit %>%
                 spread_draws(gamma[mrna_pop_num, gene_num]) %>%
                 mutate(UNIPROT = gene_map$UNIPROT[gene_num],
                        pop_mrna = mrna_pop_order[mrna_pop_num]) %>%
                 dplyr::group_by(pop_mrna, UNIPROT) %>%
                 dplyr::summarise(gamma_av = mean(gamma))
    return(gamma_info)             
}

# function that returns posterior samples of a (column normalization term)
extract_a = function(stan_fit, clusters, n_mrna){
    mrna_pop_order = paste0("Pop", 1:n_mrna)
    a_info = stan_fit %>%
             spread_draws(a_mrna[mrna_pop_num, celltype_num]) %>%
             mutate(ct = clusters[celltype_num],
                    pop_mrna = mrna_pop_order[mrna_pop_num])
    return(a_info)
}

# function that returns posterior predictive draws of mrna data set
extract_mrna_ppc = function(prep_list, stan_fit){
    mrna_ppc = stan_fit %>% # extract in vector form 
               spread_draws(mrna_sum_rep[total_obs_mrna]) %>%
               ungroup() 
               
    mrna = prep_list$mrna
    mrna_ppc = mrna %>%
               dplyr::mutate(total_obs_mrna = row_number(),
                             scaled_counts = mrna_sum/counts) %>%
               base::merge(mrna_ppc) %>% # combine posterior predictive and observed values
               mutate(mrna_bar_rep_counts = mrna_sum_rep/counts) %>% # scale total counts 
               filter(is.finite(mrna_bar_rep_counts))
    return(mrna_ppc)
}

# function to extract posterior predictive samples or protein data set
extract_protein_ppc = function(prep_list, stan_fit){
    protein_ppc = stan_fit %>%
                  spread_draws(protein_bar_rep[total_obs_protein]) %>% # pull vector of ppcs
                  ungroup() 
    protein = prep_list$protein 
    protein_ppc = protein %>% # combine with observed protein data
                  dplyr::mutate(total_obs_protein = row_number()) %>%
                  base::merge(protein_ppc)  
    return(protein_ppc)
}
