# function to run go testing
run_go_test = function(posterior_draws, organism = "human", min_size = 1, max_size = 150, sg1 = 0.01, sg2 = 0.001){
  UNIPROT_in = unique(posterior_draws$UNIPROT) # use genes present in model

  if(organism == "human"){
    # recover go groups associated with proteins fit in model
    GO_info = AnnotationDbi::select(org.Hs.eg.db, keys = unique(posterior_draws$UNIPROT), keytype = "UNIPROT", columns = "GO") %>%
      dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>%
      distinct(.keep_all = TRUE)

    # recover all proteins associated with present go groups
    GO_info2 = AnnotationDbi::select(org.Hs.eg.db, keys = unique(GO_info$GO), keytype = "GO", columns = "UNIPROT") %>%
      dplyr::group_by(GO) %>%
      dplyr::summarise(total_terms = length(unique(UNIPROT)),
                       nd = length(unique(intersect(posterior_draws$UNIPROT, UNIPROT)))) %>%
      ungroup()
  }

  if(organism == "mouse"){
    # recover go groups associated with proteins fit in model
    GO_info = AnnotationDbi::select(org.Mm.eg.db, keys = unique(posterior_draws$UNIPROT), keytype = "UNIPROT", columns = "GO") %>%
      dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>%
      distinct(.keep_all = TRUE)

    # recover all proteins associated with present go groups
    GO_info2 = AnnotationDbi::select(org.Mm.eg.db, keys = unique(GO_info$GO), keytype = "GO", columns = "UNIPROT") %>%
      dplyr::group_by(GO) %>%
      dplyr::summarise(total_terms = length(unique(UNIPROT)),
                       nd = length(unique(intersect(posterior_draws$UNIPROT, UNIPROT)))) %>%
      ungroup()
  }

  # recover go "terms" associated with groups (more detailed labels)
  GO_tab = AnnotationDbi::select(GO.db, keys = unique(GO_info$GO), keytype = "GOID", columns = "TERM") %>%
    ungroup() %>%
    dplyr::select(GOID, TERM) %>%
    distinct(.keep_all = TRUE)

  # combine previous information
  GO_info = GO_info %>%
    merge(GO_tab, by.x = "GO", by.y = "GOID") %>%
    merge(GO_info2)
  rm(list = c("GO_tab", "GO_info2"))

  # filter go groups by total number of genes and number of genes present in data
  GO_info = GO_info %>%
    filter(total_terms < max_size, nd > min_size)

  # combine go information with posterior draws for later testing
  posterior_draws_go = posterior_draws %>%
    base::merge(GO_info, by = "UNIPROT") %>%
    na.omit()

  # record size information
  size_info = posterior_draws_go %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(total_terms = unique(total_terms),
                     nd = unique(nd))

  # run test:
  test_res = posterior_draws_go %>%
    dplyr::group_by(.iteration, .chain, GO, ct) %>%
    dplyr::summarise(r_mean = mean(r, na.rm = T), # compute sample average across genes for each group, posterior draw
                     mu_mean = mean(mu, na.rm = T),
                     prot_mean = mean(prot, na.rm = T),
                     TERM = unique(TERM)) %>%
    ungroup() %>%
    dplyr::group_by(.iteration, .chain, ct) %>% # record average across genes not associated with each group
    dplyr::mutate(r_anti = (sum(r_mean, na.rm = T) - r_mean)/(n() - 1),
                  mu_anti = (sum(mu_mean, na.rm = T) - mu_mean)/(n() - 1),
                  prot_anti = (sum(prot_mean, na.rm = T) - prot_mean)/(n() - 1)) %>%
    ungroup() %>%
    dplyr::group_by(GO, ct) %>%
    dplyr::summarise(p_r = mean(r_mean < r_anti, na.rm = T), # posterior probability r
                     p_mu = mean(mu_mean < mu_anti, na.rm = T), # posterior probability mrna
                     p_prot = mean(prot_mean < prot_anti, na.rm = T), # posterior probability protein
                     TERM = unique(TERM),
                     r_off = mean(r_anti, na.rm = T), # average r across genes not in group
                     mu_off = mean(r_anti, na.rm = T), # average mu across genes not in group
                     prot_off = mean(prot_anti, na.rm = T), # average mu + r across genes not in group
                     r_av = mean(r_mean, na.rm = T), # average r in group
                     mu_av = mean(mu_mean, na.rm = T), # average mu in group
                     prot_av = mean(prot_mean, na.rm = T), # average mu + r in group
                     r_lwr = quantile(r_mean, 0.025), # 95 pct interval r
                     r_upr = quantile(r_mean, 0.975),
                     r_0 = r_lwr < 0 & r_upr > 0,
                     prt = ifelse(p_r <= (1 - p_r), p_r, (1 - p_r)), # two sided posterior probabilities
                     prm = ifelse(p_mu <= (1 - p_mu), p_mu, (1 - p_mu)),
                     prp = ifelse(p_prot <= (1 - p_prot), p_prot, (1 - p_prot))) %>%
    ungroup() %>%
    merge(size_info) %>% # bring back size information
    distinct(.keep_all = TRUE) %>%
    dplyr::group_by(ct) %>%
    arrange(prt, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
    dplyr::mutate(fdr = cummean(prt[order(prt)]), # expected proportion of false discoveries
                  sig_rt = fdr <= sg1, # expected proportion < threshold 1
                  sig_rt2 = fdr <= sg2) %>% # expected proportion < threshold 2
    arrange(prm, .by_group = T) %>%
    dplyr::mutate(fdr_m = cummean(prm[order(prm)]), # expected proportion of false discoveries for mu
                  sig_m = fdr_m <= sg1, # expected proportion < threshold 1
                  sig_m2 = fdr_m <= sg2) %>% # expected proportion < threshold 2
    arrange(prp, .by_group = T) %>%
    dplyr::mutate(fdr_p = cummean(prp[order(prp)]), # expected proportion of false discoveries for mu + R
                  sig_p = fdr_p <= sg1, # expected proportion < threshold 1
                  sig_p2 = fdr_p <= sg2) # expected proportion < threshold 2

  list(posterior_draws_go = posterior_draws_go, test_res = test_res)
}

# output list of data frames with summary of significant results
output_significant_groups = function(test_res, sg = 0.01){
  go_info = test_res %>% # select relevant summary columns
    dplyr::select(GO, ct, TERM, total_terms, nd,
                  prt, prm, prp,
                  fdr, fdr_m, fdr_p,
                  r_av, r_off,
                  mu_av, mu_off,
                  prot_av, prot_off,
                  sig_rt, sig_m, sig_p) %>%
    distinct(.keep_all = TRUE) %>%
    arrange(fdr) %>% # arrange by false discovery rate
    group_split(ct) # seprate list by cell type

  cn = c("GO", "CT", "TERM", "Total Genes", "Genes in Data",
         "Ratio p", "mRNA p", "Prot p",
         "Ratio efd","mRNA efd", "Protein efd",
         "Ratio GO", "Ratio Other",
         "mRNA GO", "mRNA Other",
         "Protein GO", "Protein Other",
         "Significant Ratio", "Significant mRNA", "Significant Protein")
  go_info_filtered = lapply(go_info, function(x) filter(x, fdr <= sg)) # filter to significant genes
  go_info_filtered = lapply(go_info_filtered, setNames, nm = cn) # set column names
  return(go_info_filtered)
}

# output table with number of significant groups associated with each parameter
go_summary_table = function(test_res){
  df = test_res %>%
    dplyr::group_by(ct) %>% # separate for each cell type
    dplyr::summarise(r_sig = sum(sig_rt), # total number of significant groups threshold 1 for R
                     r_sig2 = sum(sig_rt2), # total number of significant groups threshold 2 for R
                     mu_sig = sum(sig_m), # total number of significant groups threshold 1 for mu
                     mu_sig2 = sum(sig_m2), # total number of significant groups threshold 2 for mu
                     prot_sig = sum(sig_p), # total number of significant groups threshold 1 for mu + R
                     prot_sig2 = sum(sig_p2)) # total number of significant groups threshold 2 for mu + R
  # set column names
  colnames(df) = c("Cell Type", "R: Significant at 1 Pct FDR", "R: Significant at 0.1 Pct FDR",
                   "mRNA: Significant at 1 Pct FDR", "mRNA: Significant at 0.1 Pct FDR",
                   "Protein: Significant at 1 Pct FDR", "Protein: Significant at 0.1 Pct FDR")
  return(df)
}

# return ggplot with heatmap showing sample average of posterior means for each group
plot_go_heatmap = function(test_res, go_select = NULL){
  if(is.null(go_select)){
    go_select = test_res %>% filter(significant == TRUE) %>% pull(GO) %>% sample(20) # unless list of groups are provided, random sample
  }
  test_res = test_res %>%
    filter(GO %in% go_select) # filter to select groups

  # format df wider for posterior mean of R
  cluster_df = test_res %>%
    dplyr::select(GO, ct, r_av) %>%
    pivot_wider(names_from = ct, values_from = r_av)

  # hierarchical clustering on formatted df (GO level)
  hc = hclust(dist(dplyr::select(cluster_df, -GO)))
  go_order = data.frame(clust_order = hc$order, GO = cluster_df$GO) %>%
    mutate(GO_fct = factor(GO, levels = unique(GO[clust_order])))

  # format df wider for posterior mean of R
  cluster_df = test_res %>%
    dplyr::select(GO, ct, r_av) %>%
    pivot_wider(names_from = GO, values_from = r_av)

  # hierarchical clustering on formatted df (ct level)
  hc = hclust(dist(dplyr::select(cluster_df, -ct)))
  ct_order = data.frame(celltype_order = hc$order, ct = cluster_df$ct) %>%
    mutate(ct_fct = factor(ct, levels = unique(ct[celltype_order])))

  # merge ordering for hierarchical clustering
  test_res = test_res %>% merge(go_order, by = "GO") %>% merge(ct_order, by = "ct")
  test_res = test_res %>% arrange(GO_fct)

  # generate heatmap
  g = ggplot(data = test_res, mapping = aes(x = ct, y = forcats::fct_inorder(TERM), fill = r_av)) +
    geom_tile() +
    ylab("") +
    xlab("Cell Type") +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30), position = "right") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
                         na.value = NA, name = "")
  theme(axis.text = element_text(size = 60),
        axis.title = element_text(size = 60),
        legend.key.size = unit(1, "cm"),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', color = "slategray4"),
        panel.grid.major = element_line(color = 'slategray2'),
        panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# return point plot of posterior means
plot_go_means = function(test_res){
  # point plot of posterior means of mu and mu + R for each go group
  g = ggplot(test_res, aes(x = mu_av, y = prot_av)) +
    geom_point(aes(color = significant), size = 4) +
    facet_wrap(vars(ct), ncol = 3) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_abline(intercept = 0, slope = 1, color = "red", size = 1) +
    scale_color_manual(name = "", values = c("gray47", "purple")) +
    xlab("mRNA Posterior Mean (Among GO Groups)") +
    ylab("Protein Posterior Mean (Among GO Groups)")  +
    theme(axis.title = element_text(size = 40),
          axis.text = element_text(size = 40),
          text = element_text(size = 40),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot z transformed log fold changes of observed data for selected go group
plot_zscore_groups = function(mrna_z, protein_z, test_res,
                              go, clusters = NULL, mrna_pop_label, protein_pop_label,
                              mrna_color, protein_color){
  n_pop_m = length(mrna_z)
  n_pop_p = length(protein_z)

  # filter to preselected clusters if provided
  if(!is.null(clusters)){
    mrna_z = lapply(mrna_z, function(x) filter(x, ct %in% clusters))
    protein_z = lapply(protein_z, function(x) filter(x, ct %in% clusters))
    test_res = test_res %>% filter(ct %in% clusters)
  }

  # for each mrna data set, filter to relevant go groups
  mrna_z = lapply(mrna_z, function(x) filter(x, GO == go)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  # for each protein data set, filter to relevant go groups
  protein_z = lapply(protein_z, function(x) filter(x, GO == go)) %>%
    lapply(., function(x) arrange(x, ct)) %>%
    lapply(., function(x) mutate(x, ct_factor = factor(ct, levels = clusters)))

  # filter test results to relevant go groups
  test_res = test_res %>%
    filter(GO == go) %>%
    arrange(ct) %>%
    mutate(ct_factor = factor(ct, levels = clusters))

  # recover gene ontology group term name for relevant groups
  term_name_select = mrna_z %>% rlist::list.rbind() %>% filter(GO == go) %>% pull(TERM) %>% unique()
  if(length(term_name_select) >= 1){
    term_name_select = term_name_select %>% sample(1) # if more than one term name is associated with group, randomly select one
  }
  ylab_select = paste(term_name_select, go)

  # label mrna and protein data set
  mrna_label = paste0("mRNA ", mrna_pop_label)
  protein_label = paste0("Protein ", protein_pop_label)

  # for each data set and modality, draw line plot of log fold changes for each cell type
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
    ylab("Z Transformed Log2 GO Group, Cluster Average") +
    xlab("Cell Type") +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "", values = c(mrna_color, protein_color)) +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 30),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))

  # return plot displaying posterior interval for each cell type in group
  g2 = ggplot() +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point(data = test_res, mapping = aes(x = ct_factor, y = r_av, color = significant), size = 10) +
    geom_linerange(data = test_res, mapping = aes(x = ct_factor, ymin = r_lwr, ymax = r_upr, color = significant), linewidth = 3) +
    geom_path(data = test_res, mapping = aes(x = ct_factor, y = r_av, group = GO)) +
    ggtitle(ylab_select) +
    xlab("Cell Type") +
    ylab("Ratio 95 Pct Posterior Interval") +
    scale_y_continuous(position = "right") +
    scale_color_manual(name = "Significant rPTR", values = c("gray47", "purple")) +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 30),
          legend.key.size = unit(3, "cm"),
          legend.position = "top",
          legend.justification = "right",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g2 / g1)
}

# draw tick mark plot displaying posterior mean for genes associated with selected significant go groups
plot_ticks = function(GO_dat, gene_res,
                      n = 10, selection1 = FALSE,
                      go_type = "random", seed_id = 15,
                      go = NULL, celltype, celltype_title,
                      min_x = -4.5, max_x = 4.5){

  plot_param = "r" # use ratio as plotting parameter
  x_title = "Relative Protein to RNA Ratio" # axis label

  # randomly select groups in one of two ways
  if(selection1 == FALSE){
    GO_list = GO_dat %>% pull(GO) %>% unique() %>% sample(min(7, length(unique(GO_dat$GO)))) # randomly sample among significant genes
    GO_dat = GO_dat %>%
      dplyr::filter(GO %in% GO_list) %>%
      dplyr::mutate(TERM_factor = fct_reorder(TERM, go_mean, .desc = TRUE))
  }

  if(selection1 == TRUE){
    GO_list = GO_dat %>%
      dplyr::slice_min(prt, n = min(15, length(unique(GO_dat$GO)))) %>% # slice to most extreme significant genes
      dplyr::group_by(GO) %>%
      dplyr::summarise(ng = length(unique(UNIPROT)), # number of significant genes
                       gm = mean(param_mean)) %>% # average R among group
      dplyr::slice_max(ng, n = min(35, length(unique(GO_dat$GO)))) %>% # slice to groups with highest number of genes
      dplyr::sample_n(min(7, length(unique(GO_dat$GO))), replace = T) %>%# randomly sample
      dplyr::pull(GO)
    GO_dat = GO_dat %>%
      dplyr::filter(GO %in% GO_list) %>%
      dplyr::mutate(TERM_factor = fct_reorder(TERM, go_mean, .desc = TRUE)) # arrange by group level average
  }

  # use min and max across genes associated with all groups to auto-set axis limits
  min_check = min(GO_dat$param_mean[is.finite(GO_dat$param_mean)], na.rm = T)
  max_check = max(GO_dat$param_mean[is.finite(GO_dat$param_mean)], na.rm = T)
  lim_val = max(abs(min_check), abs(max_check))

  # draw tick mark plot
  bb =  ggplot() +
    geom_point(data = GO_dat, mapping = aes(x = param_mean, y = TERM_factor, color = param_mean, group = UNIPROT),
               size = 20, stroke = 20, shape="|", alpha = 1) +
    geom_point(data = GO_dat, mapping = aes(x = go_mean, y = TERM_factor, color = go_mean),
               size = 3, stroke = 3, shape = 4) +
    xlab(paste(x_title, celltype_title)) +
    ylab("GO Term") +
    xlim(-1*(lim_val), (lim_val)) +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red",
                          midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
                          na.value = NA, name = "") +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 20)) +
    theme(text = element_text(size = 80),
          panel.background = element_rect(fill = 'white', color = "white"),
          axis.text.y = element_text(size = 30), #margin = margin(r=0)),
          axis.text.x = element_text(size = 30),
          axis.title.y = element_text(size = 40),
          axis.title.x = element_text(size = 20),
          legend.position = "bottom",
          legend.key.size = unit(3, "cm"),
          panel.grid.major = element_line(color = 'slategray1'),
          panel.grid.minor = element_line(color = 'white'))

  # marginal table displaying posterior probability associated with each group
  tbl_df = GO_dat %>% dplyr::group_by(TERM_factor) %>% dplyr::summarise(Posterior_Probability = formatC(prt, format = "e", digits = 2)) %>% arrange(TERM_factor)
  tbl = ggplot(tbl_df, aes(y = TERM_factor, x = 1)) +
    geom_text(mapping = aes(x = 1, y = TERM_factor, label = Posterior_Probability), size = 15) +
    xlab("") +
    ylab("") +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())

  # draw title panel for table
  tbl_fill = ggplot() +
    geom_text(mapping = aes(x = 1, y = 1, label = paste("Posterior", "Probability", sep = "\n")), size = 15) +
    xlab("") +
    ylab("") +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())

  # marignal density plot across all genes
  bd <- ggplot(gene_res) +
    geom_density(aes(x = r_av)) +
    xlab("") +
    ylab("") +
    xlim(-1*(lim_val), (lim_val)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank())

  b = (((bd / bb) + plot_layout(height = c(0.35, 2.5))) | ((tbl_fill / tbl) + plot_layout(height = c(0.35, 2.5)))) + plot_layout(width = c(2.5, 0.4))

  return(b)
}
