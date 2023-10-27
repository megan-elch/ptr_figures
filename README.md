# ptr_figures
Code for figures from ptr paper.

model_functions contains functions used to run model.
go_functions contains functions used for GO testing. 
plot_functions contains functions used to draw plots and extract relevant samples, and compute prerequisite statistics for plots.

figs.Rmd generates a pdf with paper figures. 
figs_genes.Rmd displays z transformed LFCs for genes significant for at least one cell type, as well as corresponding posterior R intervals. 
figs_groups.Rmd displayed z transformed LFCs averages across genes for groups significant for at least one cell type, as well as corresponding posterior R intervals. 
