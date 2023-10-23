library(dplyr)
library(rstan)
library(tidyverse)
library(boot)
library(rlist)
library(dplyr)
library(tidybayes)
library(ggExtra)
library(ggside)
library(patchwork)
library(colorspace)
library(RColorBrewer)
library(viridis)
library(kableExtra)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)

source("scripts/model_functions.R")
source("scripts/plot_functions.R")
source("scripts/go_functions.R")

seed_id = 15
project_name = "testes"
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St")

plot_path = paste0("plots/")
output_path = paste0("model_output/")
go_path = paste0("go_results/")
load(paste0("processed_data/suff_stats_", project_name, ".RData"))

############################################################################
############################ MODEL PREP AND EDA ############################ 
############################################################################

# run model prep
protein_suff = protein_suff %>% filter(pop_protein != "Pop3")
prep_list = model_prep(mrna_suff, protein_suff, n_protein = 3, clusters = clusters)
 
g = obs_mrna_prot(prep_list)
png(file = paste0(plot_path, "mrna_prot_observed.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = obs_hist_mrna(prep_list) + xlim(0, 250)
png(file = paste0(plot_path, "mrna_observed.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = obs_hist_protein(prep_list)
png(file = paste0(plot_path, "protein_observed.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = plot_protein_obs_var(prep_list, across_celltypes = T, peptide_level = T)
png(file = paste0(plot_path, "protein_obs_var_pep.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = plot_protein_obs_var(prep_list, across_celltypes = T, peptide_level = F)
png(file = paste0(plot_path, "protein_obs_var_prot.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = plot_protein_obs_var(prep_list, across_celltypes = F)
png(file = paste0(plot_path, "protein_obs_var_ct.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

############################################################################
################################ RUN MODEL ################################# 
############################################################################

# run model, save output
fit = run_model(prep_list, project_name = project_name, n_iter = 1500, n_warmup = 1000) 
stan_fit = fit$stan_fit
gene_map = fit$gene_map
save(stan_fit, gene_map, file = paste0(output_path, "fit_p3.RData"))

# extract posterior draws
posterior_draws = extract_ratios(stan_fit, gene_map, clusters)
save(posterior_draws, file =  paste0(output_path, "posterior_draws.RData"))

############################################################################
#################### PLOT FIT/ EXTRACT RELEVANT SAMPLES #################### 
############################################################################

param_sd = posterior_draws %>% 
  dplyr::group_by(UNIPROT, .iteration, .chain) %>%
  dplyr::summarise(mu_sd = sd(mu, na.rm = T),
                   r_sd = sd(r, na.rm = T),
                   prot_sd = sd(prot, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::summarise(mu_sd_av = mean(mu_sd),
                   mu_sd_lwr = quantile(mu_sd, 0.025),
                   mu_sd_upr = quantile(mu_sd, 0.975),
                   r_sd_av = mean(r_sd),
                   r_sd_lwr = quantile(r_sd, 0.025),
                   r_sd_upr = quantile(r_sd, 0.975),
                   prot_sd_av = mean(prot_sd),
                   prot_sd_upr = quantile(prot_sd, 0.025),
                   prot_sd_upr = quantile(prot_sd, 0.975))

save(param_sd, file = paste0(output_path, "param_sd.RData"))

gene_res = test_genes_centered(posterior_draws)
save(gene_res, file = paste0(output_path, "gene_res_centered.RData"))

kappa_r_info = kappa_r_info(stan_fit = stan_fit, gene_map = gene_map, clusters = clusters, n_protein = 4)
save(kappa_r_info, file = paste0(output_path, "kappa_r_info.RData"))

gamma_r_info = gamma_r_info(stan_fit = stan_fit, gene_map = gene_map, clusters = clusters, n_mrna = 2)
save(gamma_r_info, file = paste0(output_path, "gamma_r_info.RData"))

gamma_r_info = gamma_r_info %>%
               dplyr::filter(UNIPROT %in% sample)
save(gamma_r_info, file = paste0(output_path, "gamma_r_info_filtered.RData"))

gene_res = test_genes_mean(posterior_draws = posterior_draws)
save(gene_res, file = paste0(output_path, "gene_res_mean.RData"))
print("gene test complete")

significant_intervals = compute_significant_intervals(posterior_draws = posterior_draws, gene_res = gene_res)
save(significant_intervals, file = paste0(output_path, "sig_intervals.RData"))

g = plot_significant_intervals(significant_intervals)
png(file = paste0(plot_path, "significant_intervals.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("sig intervals")

load(paste0(output_path, "protein_ppc.RData"))
ppc_var_prot = compute_protein_ppc_comparison(protein_ppc)
save(ppc_var_prot, file = paste0(output_path, "ppc_var_prot.RData"))

load(paste0(output_path, "mrna_ppc.RData"))
ppc_var_mrna = compute_mrna_ppc_comparison(mrna_ppc)
save(ppc_var_mrna, file = paste0(output_path, "ppc_var_mrna.RData"))

protein_variance_ppc = compute_protein_ppc_variance(protein_ppc)
save(protein_variance_ppc, file = paste0(output_path, "protein_ppc_var_prot.RData"))
protein_variance_ppc = compute_protein_ppc_variance(protein_ppc, protein_level = F)
save(protein_variance_ppc, file = paste0(output_path, "protein_ppc_var_pep.RData"))

fit_sig = output_significant_genes(gene_res = gene_res)
lapply(1:length(fit_sig), function(i) write.csv(fit_sig[[i]], file=paste0(plot_path, "gene_test_", unique(fit_sig[[i]]$ct), ".csv"), row.names = FALSE))
print("significant csv outputted")
 
g = plot_posterior_means(gene_res = gene_res)
png(file = paste0(plot_path, "mu_prot_points.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mu prot points")
 
g = plot_cor_across_ct(posterior_draws = posterior_draws)
png(file = paste0(plot_path, "cor_across_ct.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("correlation across cell types")

g = plot_cor_across_genes(posterior_draws = posterior_draws)
png(file = paste0(plot_path, "cor_across_genes.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("correlation across genes")

scale_means = extract_scale(stan_fit = stan_fit, gene_map = gene_map)
save(scale_means, file = paste0(output_path, "scale_means.RData"))
g = plot_scale_param(scale_means)
png(file = paste0(plot_path, "scale_means.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("scale parameter")

g = plot_mu_r_cov(posterior_draws = posterior_draws, scale_means = scale_means)
png(file = paste0(plot_path, "mu_r_cor.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mu r correlation")

g = plot_mu_r_cov(posterior_draws = posterior_draws, scale_means = scale_means, type = "cov")
png(file = paste0(plot_path, "mu_r_cov.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mu r covariance")

cor_comparison = compute_correlation_comparison(gene_res = gene_res, posterior_draws = posterior_draws, prep_list = prep_list)
save(cor_comparison, file = paste0(output_path, "cor_comparison.RData"))

g = plot_correlation_comparison_med(cor_comparison = cor_comparison)
png(file = paste0(plot_path, "cor_comparison_med.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("correlation comparison median")

g = plot_correlation_comparison_grid(cor_comparison = cor_comparison)
png(file = paste0(plot_path, "cor_comparison_grid.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("correlation comparison grid")

b_info = extract_b(stan_fit = stan_fit, clusters = clusters, n_protein = 4)
save(b_info, file = paste0(output_path, "b_info.RData"))
g = plot_b(b_info)
png(file = paste0(plot_path, "b_info.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("protein b")

kappa_info = extract_kappa(stan_fit = stan_fit, gene_map = gene_map, n_protein = 4)
save(kappa_info, file = paste0(output_path, "kappa_info.RData"))
g = plot_kappa(kappa_info, n_protein = 4)
png(file = paste0(plot_path, "kappa_info.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("kappa")
 
gamma_info = extract_gamma(stan_fit = stan_fit, gene_map = gene_map, n_mrna = 2)
save(gamma_info, file = paste0(output_path, "gamma_info.RData"))

g = plot_gamma(gamma_info, n_mrna = 2)
png(file = paste0(plot_path, "gamma_info.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("gamma")

a_info = extract_a(stan_fit = stan_fit, clusters = clusters, n_mrna = 2)
save(a_info, file = paste0(output_path, "a_info.RData"))

g = plot_a(a_info, n_mrna = 2)
png(file = paste0(plot_path, "a_info.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mrna a")

mrna_res_info = compute_mrna_residuals(prep_list = prep_list, stan_fit = stan_fit)
g = plot_mrna_residuals(mrna_res_info)
png(file = paste0(plot_path, "mrna_residuals.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mrna residuals")

protein_res_info = compute_protein_residuals(prep_list = prep_list, stan_fit = stan_fit)
g = plot_protein_residuals(protein_res_info)
png(file = paste0(plot_path, "protein_residuals.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("protein residuals")

mrna_ppc = extract_mrna_ppc(prep_list = prep_list, stan_fit = stan_fit)
save(mrna_ppc, file = paste0(output_path, "mrna_ppc.RData"))
g = plot_mrna_ppc_variance(mrna_ppc = mrna_ppc, rescale = F, viridis_type = "A")
png(file = paste0(plot_path, "mrna_ppc_var1.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mrna ppc var color scheme 1")

mrna_variance_ppc = compute_mrna_ppc_variance(mrna_ppc)
save(mrna_variance_ppc, file = paste0(output_path, "mrna_ppc_var.RData"))

g = plot_mrna_ppc_variance(mrna_ppc = mrna_ppc, rescale = F, viridis_type = "D")
png(file = paste0(plot_path, "mrna_ppc_var2.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mrna ppc var color scheme 2")

g = plot_mrna_ppc_variance(mrna_ppc = mrna_ppc, rescale = T, viridis_type = "A")
png(file = paste0(plot_path, "mrna_ppc_var3.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("mrna ppc var scaled")

protein_ppc = extract_protein_ppc(prep_list = prep_list, stan_fit = stan_fit)
save(protein_ppc, file = paste0(output_path, "protein_ppc.RData"))
g = plot_protein_ppc_variance(protein_ppc = protein_ppc, viridis_type = "A")
png(file = paste0(plot_path, "protein_ppc1.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("protein ppc var color scheme 1")

protein_variance_ppc = compute_protein_ppc_variance(protein_ppc)
save(protein_variance_ppc, file = paste0(output_path, "protein_ppc_var.RData"))

g = plot_protein_ppc_variance(protein_ppc = protein_ppc, viridis_type = "D")
png(file = paste0(plot_path, "protein_ppc2.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()
print("protein ppc var color scheme 2")

# compute and draw posterior predictive correlations 
ppc_cor_obj = compute_ppc_correlations(mrna_ppc, protein_ppc)
save(ppc_cor_obj, file = paste0(output_path, "ppc_cor.RData"))

g = draw_ppc_correlations(ppc_cor_obj)
png(file = paste0(plot_path, "ppc_cor_comparison.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = draw_ppc_correlation_density(ppc_cor_obj, mrna_df = "Pop1")
png(file = paste0(plot_path, "ppc_cor_density_m1.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

g = draw_ppc_correlation_density(ppc_cor_obj, mrna_df = "Pop2")
png(file = paste0(plot_path, "ppc_cor_density_m2.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

############################################################################
########################### GO TESTING AND PLOTS ########################### 
############################################################################

go_results = run_go_test_means(posterior_draws = posterior_draws)
save(go_results, file = paste0(go_path, "go_results_mean.RData"))
posterior_draws_go = go_results$posterior_draws_go
save(posterior_draws_go, file = paste0(go_path, "go_results_draws_mean.RData"))
test_res = go_results$test_res
save(test_res, file = paste0(go_path, "go_test_results_mean.RData"))

for(i in 1:length(clusters)){
    GO_dat = compute_ticks_draws(posterior_draws_go = posterior_draws_go, test_res = test_res, celltype = clusters[i])
    save(GO_dat, file = paste0(go_path, "GO_dat_", clusters[i], ".RData"))    
}

g = output_significant_groups(test_res = test_res, sg = 0.01)
lapply(1:length(g), function(x) write.csv(g[[x]], file=paste0(plot_path, "go_test1_", unique(g[[x]]["CT"]), ".csv"), row.names = FALSE))

g = go_summary_table(test_res = test_res)
kable(g) %>% 
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F) %>% 
  save_kable(file = paste0(plot_path, "/go_results_summary.png")) 

g = plot_go_means(test_res = test_res)
png(file = paste0(plot_path, "mu_prot_points_go.png"), type = "cairo", width = 1250, height = 750)
print(g)
dev.off()

celltype_title = c("Endothelial Cells", "Leydig Cells", "Peritubular Myoid Cells", "Spermatocyte Cells", "Spermatagonia Cells", "Spermatids")
for(i in 1:length(clusters)){
  g1 = plot_ticks(posterior_draws_go = posterior_draws_go, test_res = test_res, gene_res = gene_res, selection1 = FALSE, celltype = clusters[i], celltype_title = celltype_title[i])
  png(file = paste0(plot_path, "r_tick1_", clusters[i], ".png"), type = "cairo", width = 750, height = 750)
  print(g1)
  dev.off()
  
  g2 = plot_ticks(posterior_draws_go = posterior_draws_go, test_res = test_res, gene_res = gene_res, selection1 = TRUE, celltype = clusters[i], celltype_title = celltype_title[i])
  png(file = paste0(plot_path, "r_tick2_", clusters[i], ".png"), type = "cairo", width = 750, height = 750)
  print(g2)
  dev.off()
}


