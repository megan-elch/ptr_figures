# source model functions
source("scripts/model_functions.R")

# point plot of mrna, protein log fold changes for each cell type
obs_mrna_prot = function(prep_list){
  mrna = prep_list$mrna
  protein = prep_list$protein

  # compute protein-level averages based on peptide-level log2 intensities
  protein = protein %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(prot_av = mean(pep_av, na.rm = T))

  # compute mrna lfcs
  mrna = mrna %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(mrna_av = mean(log2(mrna_av), na.rm = T))

  df = merge(mrna, protein)

  # draw plots
  g = ggplot() +
    geom_point(data = df, mapping = aes(x = mrna_av, y = prot_av)) +
    facet_wrap(~ct, ncol = 3) +
    xlab("Mean Log2 Transcript Count") +
    ylab("Mean Log2 Protein Intensity") +
    ggtitle("mRNA and Protein Observed Averages") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# histogram of mrna averages
obs_hist_mrna = function(prep_list){
  mrna = prep_list$mrna
  # compute averages
  g_vline = mrna %>% ungroup() %>% dplyr::group_by(pop_mrna) %>% dplyr::summarise(obs_med = median(mrna_av, na.rm = T))

  # draw histogram
  g = ggplot(mrna, aes(x = mrna_av)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = obs_med)) +
    ggtitle("Observed Transcript Level Averages") +
    facet_wrap(pop_mrna ~ ct) +
    xlab("Transcript Average Across Cells") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# histgram of protein level averages
obs_hist_protein = function(prep_list){
  protein = prep_list$protein
  # compute averages
  g_vline = protein %>% ungroup() %>% dplyr::group_by(pop_protein) %>% dplyr::summarise(obs_med = median(pep_av, na.rm = T))

  # draw histogram
  g = ggplot(protein, aes(x = pep_av)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = obs_med)) +
    ggtitle("Observed Peptide Level Averages") +
    facet_wrap(pop_protein ~ ct) +
    xlab("Peptide Level Average") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot protein observed variances (across cell types or proteins)
plot_protein_obs_var = function(prep_list, across_celltypes = T, peptide_level = T){
  protein = prep_list$protein
  # compute across cell types variance if true (at peptide or protein levels)
  if(across_celltypes == TRUE){
    if(peptide_level == TRUE){
      df = protein %>%
        dplyr::group_by(pep, pop_protein) %>%
        dplyr::summarise(obs_var = var(pep_av, na.rm = T))
    }
    if(peptide_level == FALSE){
      df = protein %>%
        dplyr::group_by(UNIPROT, pop_protein) %>%
        dplyr::summarise(obs_var = var(pep_av, na.rm = T))
    }
    g_title = "Observed Variance Across Cell Types" # plot title
  }

  # compute across protein variances if false
  if(across_celltypes == FALSE){
    df = protein %>%
      dplyr::group_by(ct, pop_protein) %>%
      dplyr::summarise(obs_var = var(pep_av, na.rm = T))

    g_title = "Observed Variance Across Genes" # plot title
  }

  # vertical line at median
  g_vline = df %>% ungroup() %>% dplyr::group_by(pop_protein) %>% dplyr::summarise(var_med = median(obs_var, na.rm = T))

  # draw histogram
  g = ggplot(df, aes(x = obs_var)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = var_med)) +
    ggtitle(g_title) +
    facet_wrap(vars(pop_protein)) +
    xlab("Observed Variance") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

prepare_score = function(df, df_info, clusters = NULL, type = "mrna", protein_label = "prot", peptide_label = "pep"){
  df_info = na.omit(df_info)
  if(tolower(type) == "mrna"){
    df = df %>%
      as.matrix() %>%
      data.frame()
    df = df[,intersect(df_info$id, colnames(df))]
    df$SYMBOL = rownames(df)
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id") %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      base::merge(df_info, by = "id")
    means = means %>%
      # dplyr::group_by(SYMBOL) %>%
      # dplyr::mutate(value = value/mean(value, na.rm = T)) %>%
      dplyr::group_by(SYMBOL, ct) %>%
      dplyr::summarise(mrna_zscore = mean(log2(value + 1), na.rm = T)) %>% # mean z transformed log2 transcript count
      ungroup() %>%
      dplyr::group_by(SYMBOL) %>%
      # dplyr::mutate(mrna_zscore = scale(mrna_zscore)[,1]) # REMOVE FOR Z3
      dplyr::mutate(mrna_zscore = mrna_zscore - mean(mrna_zscore, na.rm = T)) ## PUT BACK FOR Z3 THIS LINE ONLY
    # dplyr::group_by(ct) %>%
    # dplyr::mutate(mrna_zscore = log2((mrna_zscore + 1)*sum(mrna_zscore + 1, na.rm = T))) %>%
    # ungroup()
  }
  if(tolower(type) == "protein"){
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id")
    means[["UNIPROT"]] = means[[protein_label]]
    means[["pep"]] = means[[peptide_label]]
    means = means %>%
      base::merge(df_info) %>%
      dplyr::group_by(UNIPROT, id) %>%
      dplyr::mutate(value = mean(value, na.rm = T)) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      # dplyr::group_by(UNIPROT) %>%
      # dplyr::mutate(value = value - mean(value, na.rm = T)) %>%
      # ungroup() %>%
      # dplyr::mutate(value = ifelse(is.finite(value), value, NA)) %>%
      dplyr::group_by(UNIPROT, ct) %>%
      dplyr::summarise(protein_zscore = mean(value, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT) %>%
      # dplyr::mutate(protein_zscore = scale(protein_zscore)[,1]) # REMOVE FOR Z3
      dplyr::mutate(protein_zscore = protein_zscore - mean(protein_zscore, na.rm = T)) # PUT BACK FOR Z3
  }
  return(means)
}

prepare_zscore2 = function(df, df_info, clusters = NULL, type = "mrna", protein_label = "prot", peptide_label = "pep"){
  df_info = na.omit(df_info)
  if(tolower(type) == "mrna"){
    df = df %>%
      as.matrix() %>%
      data.frame()
    df = df[,intersect(df_info$id, colnames(df))]
    df$SYMBOL = rownames(df)
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id") %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      base::merge(df_info, by = "id")

    means = means %>%
      dplyr::group_by(SYMBOL, ct) %>%
      dplyr::summarise(mrna_zscore = sum(value, na.rm = T)) %>% # mean z transformed log2 transcript count
      ungroup() %>%
      dplyr::group_by(ct) %>%
      dplyr::mutate(mrna_zscore = log2((mrna_zscore + 1)/sum(mrna_zscore + 1, na.rm = T))) %>%
      ungroup() %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::mutate(mrna_zscore = scale(mrna_zscore)[,1])
  }
  if(tolower(type) == "protein"){
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id")
    means[["UNIPROT"]] = means[[protein_label]]
    means[["pep"]] = means[[peptide_label]]
    means = means %>%
      base::merge(df_info) %>%
      dplyr::group_by(UNIPROT, id) %>%
      dplyr::mutate(value = mean(value, na.rm = T)) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(centered_value = value - mean(value, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT) %>%
      dplyr::mutate(value = centered_value) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      dplyr::mutate(value = ifelse(is.finite(value), value, NA)) %>%
      dplyr::group_by(UNIPROT, ct) %>%
      dplyr::summarise(protein_zscore = mean(value, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT) %>%
      dplyr::mutate(protein_zscore = scale(protein_zscore)[,1])
  }
  return(means)
}

prepare_zscore_data = function(df, df_info, clusters = NULL, type = "mrna", protein_label = "prot", peptide_label = "pep"){
  df_info = na.omit(df_info)
  if(tolower(type) == "mrna"){
    df = df %>%
      as.matrix() %>%
      data.frame()
    df = df[,intersect(df_info$id, colnames(df))]
    df$SYMBOL = rownames(df)
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id") %>%
      dplyr::mutate(value = log2(as.numeric(value) + 1)) %>%
      dplyr::mutate(value = ifelse(value == 0, NA, value))
    means = means %>% dplyr::group_by(id) %>%
      dplyr::mutate(value = value - sum(value, na.rm = T),
                    value2 = value - mean(value, na.rm = T)) %>%
      ungroup() %>%
      base::merge(df_info) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::mutate(value = scale(value)[,1],
                    value2 = scale(value2)[,1]) %>%
      ungroup() %>%
      dplyr::mutate(value = ifelse(is.finite(value), value, NA),
                    value2 = ifelse(is.finite(value2), value, NA)) %>%
      group_by(SYMBOL, ct) %>%
      dplyr::summarise(mrna_zscore = mean(value, na.rm = T),
                       mrna_zscore2 = mean(value2, na.rm = T)) %>% # mean z transformed log2 transcript count
      ungroup()
  }
  if(tolower(type) == "protein"){
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id")
    means[["UNIPROT"]] = means[[protein_label]]
    means[["pep"]] = means[[peptide_label]]
    means = means %>%
      base::merge(df_info) %>%
      dplyr::group_by(UNIPROT, id) %>%
      dplyr::mutate(value = mean(value, na.rm = T)) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(centered_value = value - mean(value, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT) %>%
      dplyr::mutate(orig_value = value,
                    value = scale(centered_value)[,1]) %>% # number of cells observed for each protein cell type
      ungroup() %>%
      dplyr::mutate(value = ifelse(is.finite(value), value, NA)) %>%
      dplyr::group_by(UNIPROT, ct) %>%
      dplyr::summarise(protein_zscore = mean(value, na.rm = T)) %>%
      ungroup()
  }
  return(means)
}

# plot fitted averages for mu or mu + R
plot_fit_means = function(gene_res, param = "mu"){
  gene_res = gene_res %>%
    arrange(ct)

  # select relevant parameter
  if(param == "mu"){
    gene_res = gene_res %>% mutate(param_av = mu_av)
  }
  if(param == "prot"){
    gene_res = gene_res %>% mutate(param_av = prot_av)
  }

  # point plot with path
  g = ggplot() +
    geom_point(data = gene_res, mapping = aes(x = fct_inorder(ct), y = param_av), size = 10, alpha = 0.5) +
    geom_path(data = gene_res, mapping = aes(x = fct_inorder(ct), y = param_av, group = UNIPROT), alpha = 0.25) +
    ggtitle("") +
    xlab("Cell Type") +
    ylab("Posteiror Mean") +
    scale_color_discrete(name = "") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 30),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot z transformed averages of mrna, protein and posterior intervals of significant genes
plot_zscore_means = function(mrna_z, protein_z, gene_res,
                             uni, organism = "human", clusters = NULL,
                             mrna_pop_label, protein_pop_label, mrna_color, protein_color){
  n_pop_m = length(mrna_z)
  n_pop_p = length(protein_z)

  # unify gene labels based on human or mouse organism
  if(organism == "human"){
    symbols_m = unique(unlist(lapply(mrna_z, function(x) pull(x, SYMBOL))))
    gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = symbols_m, keytype = "SYMBOL", columns = c("UNIPROT", "GENENAME"))
    mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  }
  if(organism == "mouse"){
    symbols_m = unique(unlist(lapply(mrna_z, function(x) pull(x, SYMBOL))))
    gene_rec = AnnotationDbi::select(org.Mm.eg.db, keys = symbols_m, keytype = "SYMBOL", columns = c("UNIPROT", "GENENAME"))
    mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  }

  # filter to pre-selected clusters (optional)
  if(!is.null(clusters)){
    mrna_z = lapply(mrna_z, function(x) filter(x, ct %in% clusters))
    protein_z = lapply(protein_z, function(x) filter(x, ct %in% clusters))
    gene_res = gene_res %>% filter(ct %in% clusters)
  }

  # filter to selected gene for each modality, arrange by cell type
  mrna_z = lapply(mrna_z, function(x) filter(x, UNIPROT == uni)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  protein_z = lapply(protein_z, function(x) filter(x, UNIPROT == uni)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  # same thing for fitted parameter data
  gene_res = gene_res %>%
    merge(gene_rec) %>%
    filter(UNIPROT == uni) %>%
    arrange(ct) %>%
    mutate(ct_factor = factor(ct, levels = clusters))

  # record whether observed value fits in 95 pct posterior interval for each modality
  mrna_z = lapply(mrna_z, function(x) merge(x, gene_res)) %>%
           lapply(., function(x) mutate(x, param_in = mrna_zscore >= mu_lwr & mrna_zscore <= mu_upr))

  protein_z = lapply(protein_z, function(x) merge(x, gene_res)) %>%
              lapply(., function(x) mutate(x, param_in = protein_zscore >= prot_lwr & protein_zscore <= prot_upr))

  # use full gene names for selected uniprot id
  gene_name_select = mrna_z %>% rlist::list.rbind() %>% filter(UNIPROT == uni) %>% pull(GENENAME) %>% unique()
  if(length(gene_name_select) >= 1){
    gene_name_select = gene_name_select %>% sample(1) # if there is more than one gene name present, randomly select one
  }
  ylab_select = paste(gene_name_select, uni) # axis label

  mrna_label = paste0("mRNA ", mrna_pop_label)
  protein_label = paste0("Protein ", protein_pop_label)

  # draw plot
  g1 = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = mrna_z[[1]], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical"), size = 10) +
    geom_path(data = mrna_z[[1]][!is.na(mrna_z[[1]]$mrna_zscore),], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical", group = mrna_label[1])) +
    geom_point(data = mrna_z[[2]], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical"), size = 10) +
    geom_path(data = mrna_z[[2]][!is.na(mrna_z[[2]]$mrna_zscore),], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical", group = mrna_label[2])) +
    geom_point(data = protein_z[[1]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[1]][!is.na(protein_z[[1]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[1])) +
    geom_point(data = protein_z[[2]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[2]][!is.na(protein_z[[2]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[2])) +
    geom_point(data = protein_z[[3]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[3]][!is.na(protein_z[[3]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[3])) +
    # geom_point(data = protein_z[[4]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    # geom_path(data = protein_z[[4]][!is.na(protein_z[[4]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[4])) +
    ylab("Z Log2 Cluster Level Average") +
    xlab("Cell Type") +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "", values = c(mrna_color, protein_color)) +
    # scale_shape_manual(name = "Observed in Param Interval", values = c(4, 1)) +
    theme(text = element_text(size = 50),
          # axis.text = element_text(size = 20),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  # draw posterior interval for ratio
  g2 = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = gene_res, mapping = aes(x = ct_factor, y = r_av, color = significant), size = 10) +
    geom_linerange(data = gene_res, mapping = aes(x = ct_factor, ymin = r_lwr, ymax = r_upr, color = significant), linewidth = 3) +
    geom_path(data = gene_res, mapping = aes(x = ct_factor, y = r_av, group = UNIPROT)) +
    ggtitle(ylab_select) +
    xlab("Cell Type") +
    ylab("Ratio 95 Pct Posterior Interval") +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "Significant rPTR", values = c("gray47", "purple")) +
    theme(text = element_text(size = 50),
          # axis.text = element_text(size = 30),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g2 / g1)
}

# plot z transformed averages of mrna, protein and posterior intervals of significant genes
plot_zscore_means_simple = function(mrna_z, protein_z, gene_res,
                             uni, organism = "human", clusters = NULL,
                             mrna_pop_label, protein_pop_label, mrna_color, protein_color){
  n_pop_m = length(mrna_z)
  n_pop_p = length(protein_z)

  # unify gene labels based on human or mouse organism
  if(organism == "human"){
    symbols_m = unique(unlist(lapply(mrna_z, function(x) pull(x, SYMBOL))))
    gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = symbols_m, keytype = "SYMBOL", columns = c("UNIPROT", "GENENAME"))
    mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  }
  if(organism == "mouse"){
    symbols_m = unique(unlist(lapply(mrna_z, function(x) pull(x, SYMBOL))))
    gene_rec = AnnotationDbi::select(org.Mm.eg.db, keys = symbols_m, keytype = "SYMBOL", columns = c("UNIPROT", "GENENAME"))
    mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  }

  # filter to pre-selected clusters (optional)
  if(!is.null(clusters)){
    mrna_z = lapply(mrna_z, function(x) filter(x, ct %in% clusters))
    protein_z = lapply(protein_z, function(x) filter(x, ct %in% clusters))
    gene_res = gene_res %>% filter(ct %in% clusters)
  }

  # filter to selected gene for each modality, arrange by cell type
  mrna_z = lapply(mrna_z, function(x) filter(x, UNIPROT == uni)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  protein_z = lapply(protein_z, function(x) filter(x, UNIPROT == uni)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  # same thing for fitted parameter data
  gene_res = gene_res %>%
    merge(gene_rec) %>%
    filter(UNIPROT == uni) %>%
    arrange(ct) %>%
    mutate(ct_factor = factor(ct, levels = clusters))

  # record whether observed value fits in 95 pct posterior interval for each modality
  mrna_z = lapply(mrna_z, function(x) merge(x, gene_res)) %>%
    lapply(., function(x) mutate(x, param_in = mrna_zscore >= mu_lwr & mrna_zscore <= mu_upr))

  protein_z = lapply(protein_z, function(x) merge(x, gene_res)) %>%
    lapply(., function(x) mutate(x, param_in = protein_zscore >= prot_lwr & protein_zscore <= prot_upr))

  # use full gene names for selected uniprot id
  gene_name_select = mrna_z %>% rlist::list.rbind() %>% filter(UNIPROT == uni) %>% pull(GENENAME) %>% unique()
  if(length(gene_name_select) >= 1){
    gene_name_select = gene_name_select %>% sample(1) # if there is more than one gene name present, randomly select one
  }
  ylab_select = paste(gene_name_select, uni) # axis label

  mrna_label = paste0("mRNA ", mrna_pop_label)
  protein_label = paste0("Protein ", protein_pop_label)

  # draw plot
  g1 = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = mrna_z[[1]], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical"), size = 10) +
    geom_path(data = mrna_z[[1]][!is.na(mrna_z[[1]]$mrna_zscore),], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical", group = mrna_label[1])) +
    geom_point(data = mrna_z[[2]], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical"), size = 10) +
    geom_path(data = mrna_z[[2]][!is.na(mrna_z[[2]]$mrna_zscore),], aes(x = ct_factor, y = mrna_zscore, color = "mRNA Empirical", group = mrna_label[2])) +
    geom_point(data = protein_z[[1]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[1]][!is.na(protein_z[[1]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[1])) +
    geom_point(data = protein_z[[2]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[2]][!is.na(protein_z[[2]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[2])) +
    geom_point(data = protein_z[[3]], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical"), size = 10) +
    geom_path(data = protein_z[[3]][!is.na(protein_z[[3]]$protein_zscore),], aes(x = ct_factor, y = protein_zscore, color = "Protein Empirical", group = protein_label[3])) +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "", values = c(mrna_color, protein_color)) +
    theme(text = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = "gray"),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'))
  # draw posterior interval for ratio
  g2 = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = gene_res, mapping = aes(x = ct_factor, y = r_av, color = significant), size = 10) +
    geom_linerange(data = gene_res, mapping = aes(x = ct_factor, ymin = r_lwr, ymax = r_upr, color = significant), linewidth = 3) +
    geom_path(data = gene_res, mapping = aes(x = ct_factor, y = r_av, group = UNIPROT)) +
    ggtitle(ylab_select) +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "", values = c("gray47", "purple")) +
    theme(text = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = "gray"),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'))
  return(g2 / g1)
}

# plot cell type, data set effect for each modality
plot_ct_effect = function(a_info, b_info, mrna_color, protein_color){
  g = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = a_info, aes(x = ct, y = a_av, color = "mRNA"), size = 10) + # mrna cell type effect
    geom_path(data = a_info, aes(x = ct, y = a_av, color = "mRNA", group = pop_mrna)) +
    geom_point(data = b_info, aes(x = ct, y = b_av, color = "Protein"), size = 10) + # protein cell type effect
    geom_path(data = b_info, aes(x = ct, y = b_av, color = "Protein", group = pop_protein)) +
    ylab("Cell Type Technical Effect") +
    xlab("Cell Type") +
    scale_color_manual(name = "", values = c(mrna_color, protein_color)) +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 30),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot across clusters correlation by across clusters variance
plot_across_clusters_var = function(var_df, cor_df, type = "protein"){

  # draw axis label
  if(type == "protein"){
    title_lab = "Protein Data Sets"
  }
  if(type == "mrna"){
    title_lab = "mRNA Data Sets"
  }

  # merge sources of variability
  df = merge(var_df, cor_df)

  # draw plot
  g = ggplot(df, mapping = aes(x = av_var, y = cors)) +
    geom_point() +
    ggtitle(title_lab) +
    xlab("Average Across Clusters Variance of Z Transformed Observations") +
    ylab("Average Across Clusters Within Modality Correlation") +
    theme(text = element_text(size = 25),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# compute across clusters mrna, protein correlation for observed lfcs (or correlation across data sets)
compute_across_clusters_correlations = function(mrna_z, protein_z, type = "multi_modal", mrna_symbol = TRUE, protein_symbol = FALSE, human = TRUE){
  # unify gene labels based on whether mrna is denoted as symbol and organism type
  if(mrna_symbol == TRUE & human == TRUE){
    symbol_list = unique(mrna_z$SYMBOL)
    gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = symbol_list, keytype = "SYMBOL", columns = c("UNIPROT"))
    mrna_z = mrna_z %>% merge(gene_rec)
  }
  if(mrna_symbol == TRUE & human == FALSE){
    symbol_list = unique(mrna_z$SYMBOL)
    gene_rec = AnnotationDbi::select(org.Mm.eg.db, keys = symbol_list, keytype = "SYMBOL", columns = c("UNIPROT"))
    mrna_z = mrna_z %>% merge(gene_rec)
  }
  if(protein_symbol == TRUE & human == TRUE){
    symbol_list = unique(protein_z$SYMBOL)
    gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = symbol_list, keytype = "SYMBOL", columns = c("UNIPROT"))
    protein_z = protein_z %>% merge(gene_rec)
  }
  if(protein_symbol == TRUE & human == FALSE){
    symbol_list = unique(protein_z$SYMBOL)
    gene_rec = AnnotationDbi::select(org.Mm.eg.db, keys = symbol_list, keytype = "SYMBOL", columns = c("UNIPROT"))
    protein_z = protein_z %>% merge(gene_rec)
  }

  # if multi modal, compute mrna, protein across modality correlation across cell types
  if(type == "multi_modal"){
    cor_df = mrna_z %>% merge(protein_z) %>% dplyr::group_by(UNIPROT) %>% dplyr::summarise(cors = cor(mrna_zscore, protein_zscore, use = "pairwise.complete.obs"))
  }

  # if mrna only, compute correlation across data sets across cell types
  if(type == "mrna"){
    cor_df = mrna_z %>% dplyr::mutate(mrna_df1_zscore = mrna_zscore) %>% dplyr::select(-mrna_zscore) %>% merge(protein_z) %>%
      dplyr::group_by(UNIPROT) %>% dplyr::summarise(cors = cor(mrna_zscore, mrna_df1_zscore, use = "pairwise.complete.obs"))
  }

  # if protein only, compute correlation across data sets across cell types
  if(type == "protein"){
    cor_df = mrna_z %>% dplyr::mutate(protein_df1_zscore = protein_zscore) %>% dplyr::select(-c(protein_zscore, protein_logmean)) %>% merge(protein_z) %>%
      dplyr::group_by(UNIPROT) %>% dplyr::summarise(cors = cor(protein_df1_zscore, protein_zscore, use = "pairwise.complete.obs"))
  }
  return(cor_df)
}

# plot histogram of across cell types correlation
plot_across_clusters_correlation = function(cor_df, color_select = "gray", group_ref = NULL){
  # vertical line at median
  cor_med = median(cor_df$cors)

  if(!is.null(group_ref)){
    cor_df = cor_df %>%
      merge(group_ref) %>%
      na.omit()

    cor_summary = cor_df %>%
      dplyr::group_by(var_group) %>%
      dplyr::summarise(cor_med = median(cors))

    g = ggplot(cor_df, aes(x = cors, y = var_group, fill = var_group)) +
      # geom_density_ridges(, alpha = 0.5) +
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.5) +
      # geom_vline(data = cor_summary, aes(xintercept = cor_med)) +
      colorspace::scale_fill_discrete_sequential(palette = "Plasma", name = "Protein |LFC|") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 0)) +
      xlab("Across Clusters Correlation") +
      theme(text = element_text(size = 35),
            axis.text = element_text(size = 20),
            legend.key.size = unit(1.5, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }

  if(is.null(group_ref)){
    g = ggplot(cor_df, aes(x = cors)) +
      geom_histogram(bins = 50, fill = color_select, alpha = 0.5) +
      geom_vline(xintercept = cor_med) +
      xlab("Across Clusters Correlation") +
      theme(text = element_text(size = 35),
            axis.text = element_text(size = 22),
            legend.key.size = unit(1.5, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  return(g)
}

# compute correlation across genes for all pairs of cells
compute_across_gene_correlations = function(df, df_info, cluster_select = NULL, type = "mrna"){
  # function to compute pairwise correlation
  compute_cors = function(cell_ids, df){
    cells = intersect(unique(cell_ids$id), unique(colnames(df)))
    ct = unique(cell_ids$ct)
    df = df[,cells]
    cors = cor(df, use = "pairwise.complete.obs")
    cors = cors[lower.tri(cors)]
    cor_df = data.frame(cor = cors, ct = ct)
    print(head(cor_df))
    print(mean(is.na(cor_df$cor)))
    return(cor_df)
  }

  # subset to preselected clusters optional
  if(!is.null(cluster_select)){
    df_info = df_info %>% filter(ct %in% cluster_select)
  }

  # format based on selected modality (collapse to protein level for peptides)
  if(type == "mrna"){
    df = data.frame(as.matrix(df))
  }
  if(type == "protein"){
    df = df %>%
      dplyr::select(-pep) %>%
      mutate_at(vars(-("prot")), function(x) as.numeric(as.character(x))) %>%
      dplyr::group_by(prot)  %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  }

  # apply correlation computation function
  cell_ids = df_info %>% na.omit() %>% group_split(ct)
  cor_df = lapply(cell_ids, function(x) compute_cors(x, df)) %>% list.rbind()
  return(cor_df)
}

# plot across genes correlation for supplied data frame of pairwise correlations
plot_across_gene_correlations = function(cor_df, type_label = "mRNA", color_levels, facet_col = 3, bins = 30){
  # vertical line at median
  cor_med_df = cor_df %>%
    dplyr::group_by(ct) %>%
    dplyr::summarise(cor_med = median(cor, na.rm = T))

  # draw histogram
  g = ggplot(cor_df, aes(x = cor, y = ..density.., fill = pop)) +
    geom_histogram(bins = bins, position = "identity", alpha = 0.5) +
    geom_vline(data = cor_med_df, aes(xintercept = cor_med)) +
    facet_wrap(vars(ct), ncol = facet_col, scales = "free") +
    scale_fill_manual(name = "", values = color_levels) +
    xlab(type_label) +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 19),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
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
test_genes = function(posterior_draws = NULL, stan_fit = NULL, gene_map = NULL, clusters = NULL){
  if(is.null(posterior_draws)){ # extract posterior draws if neccesary
    print("extracting ratios")
    posterior_draws = extract_ratios(stan_fit = stan_fit, gene_map = gene_map, clusters = clusters)
  }
  gene_res = posterior_draws %>%
    na.omit() %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(mu_av = mean(mu, na.rm = T), # posterior summaries mu
                     mu_lwr = quantile(mu, 0.025), na.rm = T,
                     mu_upr = quantile(mu, 0.975, na.rm = T),
                     mu_med = median(mu, na.rm = T),
                     r_av = mean(r), # posterior summaries r
                     r_lwr = quantile(r, 0.025, na.rm = T),
                     r_upr = quantile(r, 0.975, na.rm = T),
                     r_med = median(r, na.rm = T),
                     prot_av = mean(prot, na.rm = T), # posterior samples prot
                     prot_lwr = quantile(prot, 0.025, na.rm = T),
                     prot_upr = quantile(prot, 0.975, na.rm = T),
                     prot_med = median(prot, na.rm = T),
                     r_0 = r_lwr <= 0 & r_upr >= 0, # posterior interval contains 0
                     int_size = abs(r_upr - r_lwr)) %>% # interval size
    ungroup() %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::mutate(n_ct = length(intersect(unique(ct[is.finite(mu_av)]), unique(ct[is.finite(prot_av)])))) %>% # compute number of cel types observed
    ungroup() %>%
    na.omit() %>%
    mutate(significant = !r_0) # interval contains 0
  return(gene_res)
}

#  compute 95 pct interval for mrna, protein correlation for significant genes
compute_significant_intervals = function(posterior_draws, gene_res){

  # filter to significant genes
  genes = gene_res %>%
    dplyr::filter(significant == TRUE) %>%
    pull(UNIPROT) %>%
    unique()

  # compute across cell types correlation of mrna, protein for each posterior draw, genes
  sig_intervals = posterior_draws %>%
    dplyr::filter(UNIPROT %in% genes) %>%
    dplyr::group_by(UNIPROT, .iteration, .chain) %>%
    dplyr::summarise(mu_prot_cor = cor(mu, prot, use = "pairwise.complete.obs"),
                     n_ct = length(unique(ct))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(UNIPROT) %>% # for each gene, compute posterior 95 pct interval of correlation
    dplyr::summarise(cor_lwr = quantile(mu_prot_cor, 0.025),
                     cor_upr = quantile(mu_prot_cor, 0.975),
                     cor_med = median(mu_prot_cor),
                     n_ct = unique(n_ct)) %>%
    ungroup()
  return(sig_intervals)
}

# plot posterior interval for mrna, protein correlation
plot_significant_intervals = function(sig_intervals, gene_res, ct_select, organism = "human", sig_color = "purple"){
  # recover gene names based on organism
  if(organism == "human"){
    gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = unique(sig_intervals$UNIPROT), keytype = "UNIPROT", columns = c("GENENAME"))
  }
  if(organism == "mouse"){
    gene_rec = AnnotationDbi::select(org.Mm.eg.db, keys = unique(sig_intervals$UNIPROT), keytype = "UNIPROT", columns = c("GENENAME"))
  }

  # select genes based on selected cell type and significance
  ct_sig = gene_res %>% filter(significant == TRUE) %>% filter(ct == ct_select) %>% pull(UNIPROT) %>% unique()
  # obtain corresponding intervals
  sig_intervals = sig_intervals %>% filter(UNIPROT %in% ct_sig) %>% merge(gene_rec) %>% arrange(cor_med)

  # draw significant intervals, horizontal line at 0, point at median
  g = ggplot(data = sig_intervals, mapping = aes(x = fct_inorder(GENENAME), y = cor_med)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(color = "purple") +
    geom_linerange(aes(ymin = cor_lwr, ymax = cor_upr), color = sig_color) +
    # ggtitle("mRNA, Protein Across Cell Type Correlation Posterior Intervals") +
    ylab(paste(ct_select, "95 Pct Across Cell Types Correlation")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 27)) +
    xlab("") +
    # scale_color_continuous(name = "Number Of Cell Types Observed") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# function to return list of dfs with significant genes that can be outputted as csv
# input tested gene information or set of posterior draws
output_significant_genes = function(gene_res = NULL, posterior_draws = NULL){
  if(is.null(gene_res)){
    gene_res = test_genes(posterior_draws)
  }

  fit_sig = filter(gene_res, significant == TRUE) # filter to significant genes
  print(head(fit_sig))

  fit_sig = fit_sig %>%
    dplyr::group_split(ct)  # split by cell type
  ct_names = unlist(lapply(fit_sig, function(x) unique(x["ct"])))
  names(fit_sig) = ct_names
  return(fit_sig)
}

# write out as csv after:
# lapply(1:length(fit_info_sig), function(i) write.csv(fit_info_sig[[i]], file=paste0(plot_path, unique(fit_info_sig[[i]]$ct), ".csv"), row.names = FALSE))

# plot mu + r by mu, with color corresponding to significant test results
# one of "none", "mrna" or "protein" for posterior_interval (mrna or protein gives 95 pct mu and mu + r interval, respectively)
plot_posterior_means = function(gene_res, posterior_interval = "none", facet_col = 3, sig_color = "purple", test_name = "Significant Nonzero"){
  # draw point plot with reference diagonal, line at zero
  g = ggplot() +
    geom_point(data = gene_res, mapping = aes(x = mu_av, y = prot_av, color = significant, size = significant, alpha = significant)) +
    geom_abline(color = "red") +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    facet_wrap(~ct, ncol = facet_col) +
    xlab("mRNA Posterior Mean") +
    ylab("Protein Posterior Mean") +
    # ggtitle("mRNA and Protein Posterior Means by Gene and Celltype") +
    scale_color_manual(name = test_name, values = c("gray47", sig_color)) +
    scale_size_manual(values = c(3, 4.5), guide = "none") +
    scale_alpha_manual(values = c(0.45, 0.75), guide = "none") +
    theme(text = element_text(size = 35),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  # can optionally add error bars to x or y directions
  if(posterior_interval == "protein"){
    g = g + geom_errorbar(data = gene_res, aes(x = mu_av, y = prot_av, ymin = prot_lwr, ymax = prot_upr, color = significant))
  }

  if(posterior_interval == "mrna"){
    g = g + geom_errorbar(data = gene_res, aes(x = mu_av, y = prot_av, xmin = mu_lwr, xmax = mu_upr, color = significant))
  }
  return(g)
}

# plot correlation across cell types: input posterior samples of mu, r, prot
plot_cor_across_ct = function(posterior_draws, x_lab = 0.5, y_lab = 900000){
  fit_cor = posterior_draws %>%
    dplyr::group_by(UNIPROT, .iteration, .chain) %>%
    dplyr::summarise(mu_prot_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>% # across ct correlation for each draw
    ungroup()

  # compute summary information
  cor_df = fit_cor %>%
    ungroup() %>%
    dplyr::summarise(cor_med = median(mu_prot_cor, na.rm = T),# compute median for vertical line
                     cor_lab = paste("Median =", round(cor_med, 3)), x = x_lab, y = y_lab)
  # histogram of correlations
  g = ggplot(fit_cor, aes(x = mu_prot_cor)) +
    geom_histogram(bins = 100) +
    geom_vline(data = cor_df, aes(xintercept = cor_med)) +
    geom_label(data = cor_df, aes(x = x, y = y, label = cor_lab)) +
    ggtitle("mRNA, Protein Correlation Across Cell Types") +
    xlab("Across Cell Types Correlation of Posterior Samples") +
    theme(text = element_text(size = 35),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot correlation across genes: input posterior samples of mu, r, prot
plot_cor_across_genes = function(posterior_draws, y_lab = 450){
  # compute mu, prot correlation across genes for each cell type, posterior draw
  fit_cor = posterior_draws %>%
    dplyr::group_by(ct, .iteration, .chain) %>% # compute correlation across cell types for each draw
    dplyr::summarise(mu_prot_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>%
    ungroup()

  # summarise with median
  cor_df = fit_cor %>%
    ungroup() %>%
    dplyr::group_by(ct) %>%
    dplyr::summarise(cor_med = median(mu_prot_cor, na.rm = T), cor_lab = paste("Median =", round(cor_med, 3))) # compute median for vertical line

  # plot histogram
  g = ggplot(fit_cor, aes(x = mu_prot_cor)) +
    geom_histogram(bins = 100) +
    ylab("Count") +
    geom_vline(data = cor_df, aes(xintercept = cor_med)) + # vertical line
    facet_wrap(vars(ct), ncol = 3) +
    geom_label(data = cor_df, aes(x = cor_med - 0.025, y = y_lab, label = cor_lab)) +
    ggtitle("mRNA, Protein Correlation Across Genes") +
    xlab("Across Genes Correlation of Posterior Samples") +
    theme(text = element_text(size = 40),
          legend.key.size = unit(1.5, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  return(g)
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

# input posterior means for scale parameter or scale posterior means
plot_scale_param = function(scale_means = NULL, stan_fit = NULL, gene_map = NULL){
  # if scale posterior means aren't provided, compute posterior means
  if(is.null(scale_means)){
    scale_means = extract_scale(stan_fit, gene_map)
  }
  # draw plot
  g = ggplot(scale_means, aes(x = scale_av)) +
    geom_density() +
    xlab("Scale Parameter Posterior Mean") +
    ggtitle("Gene-Level Scale Parameter Posterior Mean") +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 20),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  return(g)
}

# plot covariance or correlation of mu and r, positive correlation may indicate issues with scale parameter
# input posterior draws of mu, r, and posterior means of scale parameter, select covariance or correlation
plot_mu_r_cov = function(posterior_draws = NULL, scale_means = NULL, type = "cor", stan_fit = NULL, gene_map = NULL, clusters = NULL){
  # extract ratio, sccale parameter if needed
  if(is.null(posterior_draws)){
    posterior_draws = extract_ratios(stan_fit, gene_map, clusters)
  }
  if(is.null(scale_means)){
    scale_means = extract_scale(stan_fit, gene_map)
  }
  # compute mu, r, correlation and covariance, identify posterior interval
  mu_r_info = posterior_draws %>%
    dplyr::group_by(UNIPROT, .iteration, .chain) %>%
    dplyr::summarise(mu_r_cov = cov(mu, r),
                     mu_r_cor = cor(mu, r),
                     var_prot = log(var(prot))) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(cov_med = median(mu_r_cov),
                     cov_lwr = quantile(mu_r_cov, 0.025),
                     cov_upr = quantile(mu_r_cov, 0.975),
                     cor_med = median(mu_r_cor),
                     cor_lwr = quantile(mu_r_cor, 0.025),
                     cor_upr = quantile(mu_r_cor, 0.975),
                     var_med = median(var_prot))

  # merge posterior means of scale
  mu_r_info = mu_r_info %>% merge(scale_means)

  # draw covariance or correlation
  if(type == "cov"){
    g = ggplot(data = mu_r_info, mapping = aes(x = var_med, y = cov_med)) +
      geom_point(size = 3, aes(color = scale_av)) +
      geom_errorbar(aes(ymin = cov_lwr, ymax = cov_upr, color = scale_av)) +
      scale_color_viridis(option = "A", name = "Posterior Mean of Scale Parameter") +
      ggtitle("Mu, R Covariance") +
      xlab("Log Protein Variance") +
      ylab("Posterior Interval for Mu, R Covariance") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 20),
            legend.key.size = unit(2, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  if(type == "cor"){
    g = ggplot(data = mu_r_info, mapping = aes(x = var_med, y = cor_med, color = scale_av)) +
      geom_point(size = 3, aes(color = scale_av)) +
      geom_errorbar(aes(ymin = cor_lwr, ymax = cor_upr, color = scale_av)) +
      scale_color_viridis(option = "A", name = "Posterior Mean of Scale Parameter") +
      ggtitle("Mu, R Correlation") +
      xlab("Log Protein Variance") +
      ylab("Posterior Interval for Mu, R Correlation") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 20),
            legend.key.size = unit(2, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  return(g)
}

# extract protein column normalization term. input stan fit and vector of cluster order
extract_b = function(stan_fit, clusters, n_protein = 2){
  protein_pop_order = paste0("Pop", 1:n_protein)
  # collect posterior draws for protein cell type, data set technical term
  b_info = stan_fit %>%
    spread_draws(b_protein[protein_pop_num, celltype_num]) %>%
    mutate(ct = clusters[celltype_num],
           pop_protein = protein_pop_order[protein_pop_num])
  return(b_info)
}

# plot posterior draws of column normalization term. input b draws or stan fit with necessary information
plot_b = function(b_info = NULL, stan_fit = NULL, clusters = NULL){
  if(is.null(b_info)){
    b_info = extract_b(stan_fit, clusters)
  }
  # draw density plot of data set, cell type effect
  g = ggplot(b_info, aes(x = b_protein)) +
    geom_density() +
    facet_wrap(pop_protein ~ ct) +
    xlab("Column Normalization Term") +
    ggtitle("Posterior Distribution of b (Protein)") +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 20),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# function that returns posterior means of kappa for each gene, data set
extract_kappa = function(stan_fit, gene_map, n_protein){
  protein_pop_order = paste0("Pop", 1:n_protein)
  # collect posterior draws for protein gene, data set technical effect
  kappa_info = stan_fit %>%
    spread_draws(kappa[protein_pop_num, gene_num]) %>% # extract posterior draws
    mutate(UNIPROT = gene_map$UNIPROT[gene_num],
           pop_protein = protein_pop_order[protein_pop_num]) %>%
    dplyr::group_by(pop_protein, UNIPROT) %>%
    dplyr::summarise(kappa_av = mean(kappa)) # compute posterior mean
  return(kappa_info)
}

# function to plot posterior means of kappa
plot_kappa = function(kappa_info = NULL, stan_fit = NULL, gene_map = NULL, n_protein = NULL){
  # plot protein gene, data set technical effect
  if(is.null(kappa_info)){
    kappa_info = extract_kappa(stan_fit, gene_map, n_protein)
  }
  g = ggplot(kappa_info, aes(x = kappa_av)) +
    geom_density() +
    facet_wrap(vars(pop_protein)) +
    xlab("Row Normalization Term") +
    ggtitle("Posterior Mean of Kappa (Protein)") +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 20),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  return(g)
}

# function that returns posterior means of gamma for each gene, data set
extract_gamma = function(stan_fit, gene_map, n_mrna){
  mrna_pop_order = paste0("Pop", 1:n_mrna)
  # collect posterior samples for mrna gene, data set effect
  gamma_info = stan_fit %>%
    spread_draws(gamma[mrna_pop_num, gene_num]) %>%
    mutate(UNIPROT = gene_map$UNIPROT[gene_num],
           pop_mrna = mrna_pop_order[mrna_pop_num]) %>%
    dplyr::group_by(pop_mrna, UNIPROT) %>%
    dplyr::summarise(gamma_av = mean(gamma)) # compute posterior mean
  return(gamma_info)
}

# plot posteiror means of gamma
plot_gamma = function(gamma_info = NULL, stan_fit = NULL, gene_map = NULL, n_mrna = NULL){
  if(is.null(gamma_info)){
    gamma_info = extract_gamma(stan_fit, gene_map, n_mrna)
  }
  # draw density plot for mrna data set, gene effect
  g = ggplot(gamma_info, aes(x = gamma_av)) +
    geom_density() +
    facet_wrap(vars(pop_mrna)) +
    xlab("Row Normalization Term") +
    ggtitle("Posterior Mean of Gamma (mRNA)") +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 20),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  return(g)
}

# function that returns posterior samples of a (column normalization term)
extract_a = function(stan_fit, clusters, n_mrna){
  mrna_pop_order = paste0("Pop", 1:n_mrna)
  # extract mrna data set, cell type effect
  a_info = stan_fit %>%
    spread_draws(a_mrna[mrna_pop_num, celltype_num]) %>%
    mutate(ct = clusters[celltype_num],
           pop_mrna = mrna_pop_order[mrna_pop_num])
  return(a_info)
}

# plot density of column normalization term
plot_a = function(a_info = NULL, stan_fit = NULL, clusters = NULL, n_mrna = NULL){
  if(is.null(a_info)){
    a_info = extract_a(stan_fit, clusters, n_mrna)
  }
  # plot mrna data set, cell type effect
  g = ggplot(a_info, aes(x = a_mrna)) +
    geom_density() +
    facet_wrap(pop_mrna ~ ct) +
    xlab("Column Normalization Term") +
    ggtitle("Posterior Distribution of a (mRNA)") +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 20),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  return(g)
}

# compute residuals using observed sum and mrna fitted mean for residuals vs fitted plot
compute_mrna_residuals = function(prep_list, stan_fit){
  # observed data
  mrna = prep_list$mrna %>%
    dplyr::mutate(obs_num = row_number())
  # pull model-fit mean
  mrna_res_info = stan_fit %>%
    spread_draws(mrna_mean[obs_num])
  iter_sample = sample(unique(mrna_res_info$.iteration), 100)
  mrna_res_info = mrna_res_info %>%
    filter(.iteration %in% iter_sample) %>%
    merge(mrna) %>% # merge observed data
    dplyr::group_by(UNIPROT, ct, pop_mrna) %>% # compute posterior mean
    dplyr::summarise(mrna_observed_log = mean(log(mrna_sum), na.rm = T),
                     mrna_fitted_log = mean(log(mrna_mean), na.rm = T),
                     residual = mrna_observed_log - mrna_fitted_log) # compute residual
  return(mrna_res_info)
}

# plot mrna residuals vs fitted
plot_mrna_residuals = function(mrna_res_info){
  # compute correlation of residual and fitted
  res_lab = mrna_res_info %>%
    dplyr::group_by(pop_mrna, ct) %>%
    dplyr::summarise(x = 8,
                     y = 1,
                     cor_r = cor(mrna_fitted_log, residual, use = "pairwise.complete.obs"),
                     cor_lab = paste("Cor =", round(cor_r, 3)))

  # draw plot
  g = ggplot(data = mrna_res_info, mapping = aes(x = mrna_fitted_log, y = residual)) +
    geom_point(size = 1) +
    geom_smooth(color = "blue", method = "lm") +
    facet_grid(pop_mrna ~ ct) +
    geom_label(data = res_lab, aes(x = x, y = y-1, label = cor_lab)) +
    ggtitle("mRNA Model Residuals") +
    xlab("mRNA Model Fit Log Sum Across Cells") +
    ylab(expression(paste("log", bar(Y), " - log [", mu, "a", gamma, "]"))) +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# compute residuals using observed protein log2 averages and protein fitted mean for residuals vs fitted plot
compute_protein_residuals = function(prep_list, stan_fit){
  # observed data
  protein = prep_list$protein %>%
    dplyr::mutate(obs_num = row_number())

  # pull model-fit mean
  protein_res_info = stan_fit %>%
    spread_draws(protein_mean[obs_num])
  iter_sample = sample(unique(protein_res_info$.iteration), 100)
  protein_res_info = protein_res_info %>%
    filter(.iteration %in% iter_sample) %>%
    merge(protein) %>%
    dplyr::group_by(pep, ct, pop_protein) %>%
    dplyr::summarise(protein_observed = mean(pep_av, na.rm = T), # compute observed protein level mean
                     protein_fit = mean(protein_mean, na.rm = T), # fitted protein level mean
                     residual = protein_observed - protein_fit) # compute residual

  return(protein_res_info)
}

# plot protein residuals vs fitted
plot_protein_residuals = function(protein_res_info){
  # compute correlation of fitted and residuals
  protein_res_lab = protein_res_info %>%
    dplyr::group_by(pop_protein, ct) %>%
    dplyr::summarise(x = -2,
                     y = -3.75,
                     cor_r = cor(protein_fit, residual, use = "pairwise.complete.obs"),
                     cor_lab = paste("Cor =", round(cor_r, 3)))

  # draw plot
  g = ggplot(data = protein_res_info, mapping = aes(x = protein_fit, y = residual)) +
    geom_point(size = 1) +
    geom_smooth(color = "blue", method = "lm") +
    facet_grid(pop_protein ~ ct) +
    geom_label(data = protein_res_lab, aes(x = x, y = y-1, label = cor_lab)) +
    ggtitle("Protein Model Residuals") +
    xlab("Protein Model Fit Average Across Cells") +
    ylab(expression(paste(bar(Y), " - [s(", mu, " + R) + b + ", kappa, "]"))) +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# function that computes the across cell types correlation of mrna, protein for both observed and model fit data
compute_correlation_comparison = function(gene_res, posterior_draws, prep_list){
  fit_cor = gene_res %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(param_cor = cor(mu_av, prot_av, use = "pairwise.complete.obs")) # correlation of posterior means of mu, mu + r

  # identify posterior intervals for mu and mu + r
  fit_intervals = posterior_draws %>%
    dplyr::group_by(UNIPROT, .iteration, .chain) %>%
    dplyr::summarise(fit_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>%
    dplyr::group_by(UNIPROT) %>% # posterior interval for each gene
    dplyr::summarise(param_cor_lwr = quantile(fit_cor, 0.025, na.rm = T),
                     param_cor_upr = quantile(fit_cor, 0.975, na.rm = T))
  # use protein and mrna averages for observed correlation
  protein = prep_list$protein %>%
    dplyr::group_by(UNIPROT, ct, pop_protein) %>%
    dplyr::summarise(prot_av = mean(pep_av, na.rm = T)) # compute protein level average
  mrna = prep_list$mrna %>% dplyr::select(UNIPROT, ct, pop_mrna, mrna_av)
  cor_comparison = protein %>% # merge modalities to compute correlation
    merge(mrna) %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::group_by(UNIPROT, pop_protein) %>%
    dplyr::summarise(n_ct_pop1 = length(unique(ct[pop_mrna == "Pop1"])), # record number of cell types
                     n_ct_pop2 = length(unique(ct[pop_mrna == "Pop2"])),
                     pop1_cor = ifelse(n_ct_pop1 > 2, # observed data correlation
                                       cor(mrna_av[pop_mrna == "Pop1"], prot_av[pop_mrna == "Pop1"], use = "pairwise.complete.obs"),
                                       NA), # mrna data set 1 correlation (between modalities, across data sets)
                     pop2_cor = ifelse(n_ct_pop2 > 2, # mrna data set 2 correlation (between modalities, across data sets )
                                       cor(mrna_av[pop_mrna == "Pop2"], prot_av[pop_mrna == "Pop2"], use = "pairwise.complete.obs"),
                                       NA)) %>%
    ungroup() %>%
    na.omit() %>%
    merge(fit_cor) %>%
    merge(fit_intervals) # combine observed and model fit correlations into one df
  return(cor_comparison)
}

# returns point plot displaying observed and model-fit correlations
plot_correlation_comparison_med = function(cor_comparison){
  cor_info = cor_comparison %>% dplyr::group_by(UNIPROT) %>% dplyr::summarise(n_ct_pop1 = round(median(n_ct_pop1, na.rm = T)), # use across-data set medians
                                                                              n_ct_pop2 = round(median(n_ct_pop2, na.rm = T)),
                                                                              param_cor = median(param_cor, na.rm = T),
                                                                              pop1_cor = median(pop1_cor, na.m = T),
                                                                              pop2_cor = median(pop2_cor, na.rm = T),
                                                                              param_cor_lwr = median(param_cor_lwr, na.rm = T),
                                                                              param_cor_upr = median(param_cor_upr, na.rm = T))

  lab_df2 = cor_info %>% # compare intervals above and below y=x
    mutate(above1 = param_cor >= pop1_cor, # record whether fitted correlation is higher than observed for each data set
           above2 = param_cor >= pop2_cor)

  lab1 = lab_df2 %>%
    dplyr::mutate(yex = param_cor_lwr <= pop1_cor & param_cor_upr >= pop1_cor) %>%
    dplyr::group_by(above1) %>%
    dplyr::summarise(cover1 = mean(yex, na.rm = T),
                     l1 = ifelse(above1 == TRUE, paste("Above mRNA Pop1 Coverage =", round(cover1, 2)), paste("Below mRNA Pop1 Coverage =", round(cover1, 2))),
                     x = ifelse(above1 == TRUE, -0.5, -0.25),
                     y = ifelse(above1 == TRUE, 0, -0.75)) %>%
    ungroup()

  lab2 = lab_df2 %>%
    dplyr::mutate(yex = param_cor_lwr <= pop2_cor & param_cor_upr >= pop2_cor) %>%
    dplyr::group_by(above2) %>%
    dplyr::summarise(cover2 = mean(yex, na.rm = T),
                     l2 = ifelse(above2 == TRUE, paste("Above mRNA Pop2 Coverage =", round(cover2, 2)), paste("Below mRNA Pop2 Coverage =", round(cover2, 2))),
                     x = ifelse(above2 == TRUE, -0.5, -0.25),
                     y = ifelse(above2 == TRUE, -0.1, -0.85)) %>%
    ungroup()

  # draw main plot
  g1 = ggplot() +
    geom_point(data = cor_info,
               aes(y = param_cor, x = pop1_cor, color = n_ct_pop1), size = 3) +
    geom_point(data = cor_info,
               aes(y = param_cor, x = pop2_cor, color = n_ct_pop1), size = 3) +
    geom_abline(color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    xlab("Median Empirical Correlation") +
    ylab("Posterior Mean Correlation") +
    scale_x_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    scale_color_continuous_sequential(palette = "Viridis", name = "Number of Cell Types Observed") +
    # geom_label(data = lab1, aes(x = x, y = y, label = l1), color = "black", alpha = 1, size = 6) +
    # geom_label(data = lab2, aes(x = x, y = y, label = l2), color = "black", alpha = 1, size = 6) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text = element_text(size = 60),
          axis.title = element_text(size = 60),
          legend.key.size = unit(3, "cm"),
          legend.position = "bottom",
          legend.title = element_text(size=60),
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  # prepare data for marginal plot
  cor_info_piv = pivot_longer(cor_comparison, cols = c("pop1_cor", "pop2_cor"), values_to = "pop_cor") %>%
                 pivot_longer(cols = c("n_ct_pop1", "n_ct_pop2"), values_to = "n_ct", names_to = "n_ct_pop")
  g2_lab = cor_info_piv %>% dplyr::summarise(pop_med = median(pop_cor, na.rm = T),
                                             lab = paste0("Median = ", round(pop_med, 2)),
                                             x = pop_med,
                                             y = 0.1)
  param_med = median(cor_info$param_cor)
  g3_lab = data.frame(y = 0.5, x = 3, lab = paste0("Median = ", round(param_med, 2)))

  # draw marginal density plots
  g2 = ggplot() +
    geom_density(data = cor_info_piv, mapping = aes(x = pop_cor, group = n_ct, color = n_ct)) +
    scale_color_continuous_sequential(palette = "Viridis", name = "") +
    geom_vline(data = g2_lab, aes(xintercept = pop_med)) +
    # geom_label(data = g2_lab, aes(x = x, y = y, label = lab), color = "black", alpha = 1, size = 15) +
    scale_x_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

  # marginal plot for fitted parameter
  g3 = ggplot() +
    geom_density(data = cor_info, mapping = aes(y = param_cor)) +
    geom_hline(yintercept = param_med) +
    # geom_label(data = g3_lab, aes(x = x, y = y, label = lab), color = "black", alpha = 1, size = 15) +
    scale_y_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

  # use patchwork to combine into one plot
  g = (((g2 / g1) + plot_layout(height = c(0.35, 2.5))) | ((plot_spacer() / g3) + plot_layout(height = c(0.35, 2.5)))) + plot_layout(width = c(2.5, 0.2))
  return(g)
}

# return similar plot as previous function, but a gridded by data set instead of the median
plot_correlation_comparison_grid = function(cor_comparison){
  lab_df = cor_comparison %>%
    dplyr::group_by(pop_protein) %>% # summary information for labels
    dplyr::summarise(pop_cor_med = paste0("Median Corr Pop1, Pop2: ",
                                          round(median(pop1_cor, na.rm = T), 2),
                                          ", ",
                                          round(median(pop2_cor, na.rm = T), 2)),
                     param_cor_med = paste("Median Corr Signal:",
                                           round(median(param_cor), 2)),
                     above_cov_pop1 = paste("mRNA Pop1:", round(mean(param_cor >= pop1_cor, na.rm = T), 2)),
                     above_cov_pop2 = paste("mRNA Pop2:", round(mean(param_cor >= pop2_cor, na.rm = T), 2)),
                     below_cov_pop1 = paste("mRNA Pop1:", round(mean(param_cor <= pop1_cor, na.rm = T), 2)),
                     below_cov_pop2 = paste("mRNA Pop2:", round(mean(param_cor <= pop2_cor, na.rm = T), 2)),
                     x_pop = 0.5, y_pop = -0.55,
                     x_param = 0.5, y_param = -0.8,
                     x_a1 = -0.75, y_a1 = -0.3,
                     x_a2 = -0.75, y_a2 = -0.55,
                     x_b1 = -0.5, y_b1 = -0.75,
                     x_b2 = -0.5, y_b2 = -0.9) %>%
    ungroup()
  # draw plot
  g = ggplot() +
    geom_point(data = cor_comparison,
               aes(y = param_cor, x = pop1_cor, color = n_ct_pop1), size = 1) +
    geom_point(data = cor_comparison,
               aes(y = param_cor, x = pop2_cor, color = n_ct_pop1), size = 1) +
    geom_abline(color = "red") +
    facet_wrap(~pop_protein, ncol = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    ggtitle("Correlation Between mRNA and Protein Across Cell Types") +
    xlab("Emprirical Correlation") +
    ylab("Inferred Correlation") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    scale_color_continuous_sequential(palette = "Viridis", name = "Number of Cell Types Observed") +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# function that returns posterior predictive draws of mrna data set
extract_mrna_ppc = function(prep_list, stan_fit){
  mrna_ppc = stan_fit %>% # extract in vector form
    spread_draws(mrna_sum_rep[total_obs_mrna]) %>%
    ungroup()

  mrna = prep_list$mrna
  mrna_ppc = mrna %>%
    dplyr::mutate(total_obs_mrna = row_number(),
                  scaled_counts = mrna_sum/counts) %>% # scale mrna counts
    base::merge(mrna_ppc) %>% # combine posterior predictive and observed values
    mutate(mrna_bar_rep_counts = mrna_sum_rep/counts) %>% # scale total counts
    filter(is.finite(mrna_bar_rep_counts))
  return(mrna_ppc)
}

compute_mrna_ppc_variance = function(mrna_ppc){
  ppc_var_plot = mrna_ppc %>%
    dplyr::group_by(UNIPROT, pop_mrna, .iteration, .chain) %>%
    dplyr::summarise(ppc_var = var(mrna_bar_rep_counts), # use variance based on counts scaled by total transcript counts
                     ppc_var_rescale = var(mrna_bar_rep_counts/mean(mrna_bar_rep_counts)), # optionally rescale by mean so variances aren't so small
                     n_c = sum(n_cells),
                     av_var = var(scaled_counts), # compute variance
                     av_var_rescale = var((scaled_counts)/mean(scaled_counts))) %>% # scale variance by mean
    ungroup() %>%
    dplyr::group_by(UNIPROT, pop_mrna) %>% # record posterior predictive intervals
    dplyr::summarise(ppc_var_lwr = quantile(ppc_var, 0.025, na.rm = T),
                     ppc_var_med = median(ppc_var, na.rm = T),
                     ppc_var_upr = quantile(ppc_var, 0.975, na.rm = T),
                     av_var = median(av_var, na.rm = T),
                     ppc_var_lwr_rescale = quantile(ppc_var_rescale, 0.025, na.rm = T),
                     ppc_var_med_rescale = median(ppc_var_rescale, na.rm = T),
                     ppc_var_upr_rescale = quantile(ppc_var_rescale, 0.975, na.rm = T),
                     av_var_rescale = median(av_var_rescale, na.rm = T),
                     n_cells = unique(n_c))
  return(ppc_var_plot)
}

# plot variance across data sets for set of genes of interest
plot_protein_across_pop_var = function(ppc_var_plot, suspcious_genes){
  # summary label info
  ppc_var_plot = ppc_var_plot %>%
    dplyr::mutate(interest = UNIPROT %in% suspicious_genes)
  ppc_coverage = ppc_var_plot %>%
    dplyr::group_by(interest, ct) %>% # compute summary information for each cell type, based on whether gene is in set
    dplyr::summarise(cover = mean(av_var >= ppc_var_lwr & av_var <= ppc_var_upr, na.rm = T),
                     x = quantile(log(av_var, 0.85), na.rm = T),
                     y = quantile(log(ppc_var_med, 0.85), na.rm = T),
                     cover_lab = paste("Coverage of Obs: ", round(cover, 2)))

  # draw plot
  g = ggplot(data = ppc_var_plot, aes(x = log(av_var), y = log(ppc_var_med), color = interest)) +
    facet_wrap(vars(ct), ncol = 3) +
    geom_linerange(aes(ymin = log(ppc_var_lwr), ymax = log(ppc_var_upr), linewidth = interest)) +
    geom_abline(color = "red", alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    ggtitle("") +
    xlab("Log Observed Variance") +
    ylab("Log PP Variance") +
    scale_color_manual(values = c("lightgray", "purple")) +
    scale_linewidth_manual(values = c(0.2, 0.8)) +
    scale_alpha_manual(values = c(0.1, 1)) +
    # geom_label(data = ppc_coverage, aes(x = x, y = y, label = cover_lab, color = interest), alpha = 1, size = 8) +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 30),
          legend.key.size = unit(1, "cm"),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot variance across data sets for mrna using set of genes of interest
plot_mrna_across_pop_var = function(ppc_var_plot, rescale = F, suspicious_genes){
  # draw plot. if using rescale (by mean), plot accordingly
  if(rescale == FALSE){
    ppc_var_plot = ppc_var_plot %>%
                   dplyr::mutate(interest = UNIPROT %in% suspicious_genes) %>%
                   filter(interest == T)
    ppc_coverage = ppc_var_plot %>%
      dplyr::group_by(interest, ct) %>%
      dplyr::summarise(cover = mean(av_var >= ppc_var_lwr & av_var <= ppc_var_upr, na.rm = T),
                       x = quantile(log(av_var, 0.85), na.rm = T),
                       y = quantile(log(ppc_var_med, 0.85), na.rm = T),
                       cover_lab = paste("Coverage of Obs: ", round(cover, 2)))

    # draw plot
    g = ggplot(data = ppc_var_plot, aes(x = log(av_var), y = log(ppc_var_med), color = interest)) +
      facet_wrap(vars(ct), ncol = 3) +
      geom_linerange(aes(ymin = log(ppc_var_lwr), ymax = log(ppc_var_upr), color = interest, linewidth = interest, alpha = interest)) +
      geom_abline(color = "red", alpha = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
      scale_color_manual(values = c("lightgray", "purple")) +
      scale_alpha_manual(values = c(0.1, 1)) +
      scale_linewidth_manual(values = c(0.2, 1.8)) +
      # geom_label(data = ppc_coverage, aes(x = x, y = y, label = cover_lab, color = interest), alpha = 1, size = 8) +
      xlab("Log Observed Variance") +
      ylab("Log PP Variance") +
      scale_color_discrete(name = "Genes of Interest") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1, "cm"),
            legend.position = "none",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  # plot for variance rescaled by mean
  if(rescale == TRUE){
    g = ggplot(data = ppc_var_plot, aes(x = log(av_var_rescale), y = log(ppc_var_med_rescale), color = interest)) +
      geom_point(size = 1) +
      geom_rug() +
      geom_errorbar(aes(ymin = log(ppc_var_lwr_rescale), ymax = log(ppc_var_upr_rescale), color = interest)) +
      geom_abline(color = "red", alpha = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
      scale_color_discrete(name = "Genes of Interest") +
      xlab("Log Observed FC Variance") +
      ylab("Log PPC FC Variance") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  return(g)
}

# plot comparison of posterior predictive and observed mrna variances
plot_mrna_ppc_variance = function(ppc_var_plot, rescale = F, viridis_type = "Viridis"){
  # draw plot. if using rescale (by mean), plot accordingly
  if(rescale == FALSE){
    ppc_coverage = ppc_var_plot %>%
      dplyr::group_by(pop_mrna) %>% # compute coverage with observed and posterior predictive variances
      dplyr::summarise(cover = mean(av_var >= ppc_var_lwr & av_var <= ppc_var_upr, na.rm = T),
                       corr = cor(av_var, ppc_var_med, use = "pairwise.complete.obs")) %>%
      ungroup() %>%
      mutate(x = -30, y = -15,
             cover_lab = paste("Coverage of Obs: ", round(cover, 2)),
             corr_lab = paste("Correlation: ", round(corr, 2)))

    # draw plot
    g = ggplot(data = ppc_var_plot, aes(x = log(av_var), y = log(ppc_var_med), color = log(n_cells))) +
      geom_point(size = 1) +
      geom_rug() +
      geom_errorbar(aes(ymin = log(ppc_var_lwr), ymax = log(ppc_var_upr), color = log(n_cells))) +
      geom_abline(color = "red", alpha = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
      facet_grid(vars(pop_mrna)) +
      # geom_label(data = ppc_coverage, aes(x = x, y = y, label = cover_lab), color = "black", alpha = 1, size = 8) +
      # geom_label(data = ppc_coverage, aes(x = x, y = y-3, label = corr_lab), color = "black", alpha = 1, size = 8) +
      xlab("Log Observed Variance") +
      ylab("Log PP Variance") +
      scale_color_continuous_sequential(palette = viridis_type, name = "Log Number of Cells") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  if(rescale == TRUE){
    # plot for rescaled variance
    g = ggplot(data = ppc_var_plot, aes(x = log(av_var_rescale), y = log(ppc_var_med_rescale), color = log(n_cells))) +
      geom_point(size = 1) +
      geom_rug() +
      geom_errorbar(aes(ymin = log(ppc_var_lwr_rescale), ymax = log(ppc_var_upr_rescale), color = log(n_cells))) +
      geom_abline(color = "red", alpha = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
      facet_grid(vars(pop_mrna)) +
      xlab("Log Observed FC Variance") +
      ylab("Log PPC FC Variance") +
      scale_color_continuous_sequential(palette = viridis_type, name = "Log Number of Cells", limits = c(0, 8), oob = scales::squish) +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  return(g)
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

compute_protein_ppc_variance = function(protein_ppc, protein_level = T){
  ppc_var_plot = protein_ppc %>%
    dplyr::group_by(pep, pop_protein, .iteration, .chain) %>%
    dplyr::summarise(ppc_var = var(protein_bar_rep), # compute posterior predictive variance for each peptide, data set, draw
                     n_c = sum(n_cells),
                     var_protein = var(pep_av, na.rm = T),
                     UNIPROT = unique(UNIPROT),
                     n_c_prot = sum(n_cells_prot)) %>%
    ungroup()

  if(protein_level == TRUE){
    ppc_var_plot = ppc_var_plot %>%
      dplyr::group_by(UNIPROT, pop_protein) %>% # compute 95 pct inverals
      dplyr::summarise(ppc_var_lwr = quantile(ppc_var, 0.025, na.rm = T),
                       ppc_var_med = median(ppc_var, na.rm = T),
                       ppc_var_upr = quantile(ppc_var, 0.975, na.rm = T),
                       n_c = median(n_c_prot),
                       var_protein = mean(var_protein))

  }
  if(protein_level == FALSE){
    ppc_var_plot = ppc_var_plot %>%
      dplyr::group_by(pep, pop_protein) %>% # compute 95 pct inverals
      dplyr::summarise(ppc_var_lwr = quantile(ppc_var, 0.025, na.rm = T),
                       ppc_var_med = median(ppc_var, na.rm = T),
                       ppc_var_upr = quantile(ppc_var, 0.975, na.rm = T),
                       n_c = sum(n_c),
                       var_protein = mean(var_protein))

  }
  return(ppc_var_plot)
}

# function to plot posterior predictive protein variance. no rescale option, as this represents peptide-averages already
plot_protein_ppc_variance = function(ppc_var_plot, viridis_type = "Viridis", protein_level = T){
  # summary label info
  ppc_coverage = ppc_var_plot %>%
    dplyr::group_by(pop_protein) %>% # summarise information
    dplyr::summarise(cover = mean(var_protein >= ppc_var_lwr & var_protein <= ppc_var_upr, na.rm = T),
                     corr = cor(var_protein, ppc_var_med, use = "pairwise.complete.obs")) %>%
    ungroup() %>%
    mutate(x = -9, y = -1,
           cover_lab = paste("Coverage of Obs: ", round(cover, 2)),
           corr_lab = paste("Correlation: ", round(corr, 2)))

  ppc_var_plot = ppc_var_plot %>% mutate(n_cells = n_c)

  # draw plot
  g = ggplot(data = ppc_var_plot, aes(x = log(var_protein), y = log(ppc_var_med), color = log(n_cells))) +
    geom_point(size = 1, aes(alpha = log(n_cells))) +
    geom_rug() +
    geom_errorbar(aes(ymin = log(ppc_var_lwr), ymax = log(ppc_var_upr), color = log(n_cells), alpha = exp(log(n_c))/(1+exp(log(n_c))))) +
    scale_color_continuous_sequential(palette = viridis_type, name = "Log Number of Cells", limits = c(0, 8), oob = scales::squish) +
    scale_alpha(guide = 'none') +
    geom_abline(color = "red", alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    facet_grid(vars(pop_protein)) +
    ggtitle("") +
    xlab("Log Observed Variance") +
    ylab("Log PP Variance") +
    # geom_label(data = ppc_coverage, aes(x = x, y = y, label = cover_lab), color = "black", alpha = 1, size = 4) +
    # geom_label(data = ppc_coverage, aes(x = x, y = y-2.5, label = corr_lab), color = "black", alpha = 1, size = 4) +
    theme(text = element_text(size = 30),
          axis.text = element_text(size = 30),
          legend.key.size = unit(1, "cm"),
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# function to compute observed and posterior predictive correlations
compute_ppc_correlations = function(mrna_ppc, protein_ppc){
  # compute protein-level average of protein posterior predictive average
  protein_ppc_av = protein_ppc %>%
    dplyr::group_by(UNIPROT, ct, pop_protein, .iteration, .chain) %>%
    dplyr::summarise(protein_bar_rep = mean(protein_bar_rep, na.rm = T))

  # compute posterior predictive correlations
  ppc_cor = mrna_ppc %>%
    merge(protein_ppc_av, by = c("UNIPROT", "ct", ".iteration", ".chain")) %>%
    dplyr::group_by(UNIPROT, pop_protein, .iteration, .chain) %>%
    dplyr::summarise(pop1_ppc = ifelse(sum(pop_mrna == "Pop1") > 3,
                                       cor(mrna_bar_rep_counts[pop_mrna == "Pop1"], protein_bar_rep[pop_mrna == "Pop1"],
                                           use = "pairwise.complete.obs"),
                                       NA),
                     pop2_ppc = ifelse(sum(pop_mrna == "Pop2") > 3,
                                       cor(mrna_bar_rep_counts[pop_mrna == "Pop2"], protein_bar_rep[pop_mrna == "Pop2"],
                                           use = "pairwise.complete.obs"),
                                       NA)) %>%
    na.omit() %>%
    ungroup()

  # compute protein-level observed averages
  obs_cor = protein_ppc %>%
    dplyr::group_by(UNIPROT, pop_protein, ct) %>%
    dplyr::summarise(prot_av = mean(pep_av, na.rm = T)) %>%
    na.omit() %>%
    ungroup()

  # mrna observed averages
  mrna_obs_cor = mrna_ppc %>%
    dplyr::group_by(UNIPROT, pop_mrna, ct) %>%
    dplyr::summarise(mrna_av = mean(scaled_counts, na.rm = T)) %>%
    na.omit() %>%
    ungroup()

  # compute observed correlations
  obs_cor = obs_cor %>%
    merge(mrna_obs_cor, by = c("UNIPROT", "ct")) %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::group_by(pop_protein, UNIPROT) %>%
    dplyr::summarise(pop1_cor = ifelse(sum(pop_mrna == "Pop1") > 3,
                                       cor(mrna_av[pop_mrna == "Pop1"], prot_av[pop_mrna == "Pop1"], use = "pairwise.complete.obs"),
                                       NA),
                     pop2_cor = ifelse(sum(pop_mrna == "Pop2") > 3,
                                       cor(mrna_av[pop_mrna == "Pop2"], prot_av[pop_mrna == "Pop2"], use = "pairwise.complete.obs"),
                                       NA)) %>%
    ungroup() %>%
    na.omit()

  list(obs_cor = obs_cor, ppc_cor = ppc_cor)
}

# errorbar plot displaying 95 pct intervals for posterior predictive and observed correlations
draw_ppc_correlations = function(ppc_cor_obj){
  # extract relevant correlations
  ppc_cor = ppc_cor_obj$ppc_cor
  obs_cor = ppc_cor_obj$obs_cor

  # 95 pct intervals for posterior predictive correlations and merge observed correlations
  cor_comparison = ppc_cor %>%
    dplyr::group_by(pop_protein, UNIPROT) %>%
    dplyr::summarise(pop1_lwr = quantile(pop1_ppc, 0.025),
                     pop1_upr = quantile(pop1_ppc, 0.975),
                     pop2_lwr = quantile(pop2_ppc, 0.025),
                     pop2_upr = quantile(pop2_ppc, 0.975),
                     pop1_med = median(pop1_ppc, na.rm = T),
                     pop2_med = median(pop2_ppc, na.rm = T)) %>%
    merge(obs_cor)

  # label with coverage and correlation of medians
  g_lab = cor_comparison %>%
    dplyr::group_by(pop_protein) %>%
    dplyr::summarise(pop1_cover = (mean(pop1_cor >= pop1_lwr & pop1_cor <= pop1_upr) + mean(pop2_cor >= pop2_lwr & pop2_cor <= pop2_upr))/2,
                     pop1_coverage = paste("Av Coverage:", round(pop1_cover, 3)),
                     pop1_cor_m = (cor(pop1_med, pop1_cor, use = "pairwise.complete.obs") + cor(pop2_med, pop2_cor, use = "pairwise.complete.obs"))/2,
                     pop1_cor_lab = paste("Av Correlation:", round(pop1_cor_m, 3)))


  # draw errorbar plot of posterior predictive correlations and observed correlations
  g = ggplot(data = cor_comparison, aes(x = pop1_cor, y = pop1_med, color = "mRNA Pop1")) +
    geom_point(size = 2.5) +
    geom_rug() +
    geom_errorbar(aes(ymin = pop1_lwr, ymax = pop1_upr, color = "mRNA Pop1"), alpha = 0.5) +
    geom_point(aes(x = pop2_cor, y = pop2_med, color = "mRNA Pop2"), size = 2.5) +
    geom_rug() +
    geom_errorbar(data = cor_comparison, aes(ymin = pop2_lwr, ymax = pop2_upr, color = "mRNA Pop2"), alpha = 0.5) +
    scale_color_discrete(name = " ") +
    geom_abline(color = "red", alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    facet_grid(vars(pop_protein)) +
    # ggtitle("Protein Observed and Posterior Predictive Correlation (Across Celltypes)") +
    xlab("Observed Correlation Across Celltypes") +
    ylab("PP Correlation Across Celltypes") +
    geom_label(data = g_lab, aes(x = 0.8, y = -0.5, label = pop1_coverage),
               color = "black", alpha = 1, size = 4) +
    geom_label(data = g_lab, aes(x = 0.8, y = -0.8, label = pop1_cor_lab),
               color = "black", alpha = 1, size = 4) +
    theme(text = element_text(size = 40),
          axis.text = element_text(size = 40),
          legend.key.size = unit(1, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# density plot showing posterior predictive and observed correlation density plots
draw_ppc_correlation_density = function(ppc_cor_obj, mrna_df = "both"){
  # extract relevant correlations
  ppc_cor = ppc_cor_obj$ppc_cor
  obs_cor = ppc_cor_obj$obs_cor

  if(mrna_df == "both"){
    ppc_cor = ppc_cor %>%
              dplyr::mutate(av_ppc = ifelse(!is.na(pop1_ppc + pop2_ppc), (pop1_ppc + pop2_ppc)/2,
                                            ifelse(is.na(pop1_ppc), pop2_ppc, pop1_ppc)))
    obs_cor = obs_cor %>%
              dplyr::mutate(av_obs = ifelse(!is.na(pop1_cor + pop2_cor), (pop1_cor + pop2_cor)/2,
                                            ifelse(is.na(pop1_cor), pop2_cor, pop1_cor)))
    g = ggplot() +
      geom_density(data = ppc_cor, mapping = aes(x = av_ppc, color = "PP", group = .iteration, y = ..density..), position = "dodge", alpha = 0.1, size = 0.5) +
      geom_density(data = obs_cor, mapping = aes(x = av_obs, color = "Observed", y = ..density..), position = "dodge", size = 1) +
      scale_color_brewer(palette = "Paired", limits = c("PP", "Observed"), name = "") +
      facet_wrap(~pop_protein) +
      xlab("Across Cell Types Correlation") +
      ylab("Density") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1.5, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'white'))
  }

  # density plot for mRNA dataset 1
  if(mrna_df == "Pop1"){
    g = ggplot() +
      geom_density(data = ppc_cor, mapping = aes(x = pop1_ppc, color = "Pop1 PP", group = .iteration, y = ..density..), position = "dodge", alpha = 0.1, size = 0.5) +
      geom_density(data = obs_cor, mapping = aes(x = pop1_cor, color = "Pop1 Observed", y = ..density..), position = "dodge", size = 1) +
      scale_color_brewer(palette = "Paired", limits = c("Pop1 PP", "Pop1 Observed"), name = "") +
      facet_wrap(~pop_protein) +
     # ggtitle(paste0("Posterior Predictive and Observed mRNA ", mrna_df, ", Protein Correlation Across Cell Types")) +
      xlab("Across Cell Types Correlation") +
      ylab("Density") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1.5, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'white'))
  }

  if(mrna_df == "Pop2"){
    g = ggplot() +
      geom_density(data = ppc_cor, mapping = aes(x = pop2_ppc, color = "Pop2 PP", group = .iteration, y = ..density..), position = "dodge", alpha = 0.1, size = 0.5) +
      geom_density(data = obs_cor, mapping = aes(x = pop2_cor, color = "Pop2 Observed", y = ..density..), position = "dodge", size = 1) +
      scale_color_brewer(palette = "Paired", limits = c("Pop2 PP", "Pop2 Observed"), name = "") +
      facet_wrap(~pop_protein) +
      # ggtitle(paste0("Posterior Predictive and Observed mRNA ", mrna_df, ", Protein Correlation Across Cell Types")) +
      xlab("Across Cell Types Correlation") +
      ylab("Density") +
      theme(text = element_text(size = 30),
            axis.text = element_text(size = 30),
            legend.key.size = unit(1.5, "cm"),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'white'))
  }
  return(g)
}
