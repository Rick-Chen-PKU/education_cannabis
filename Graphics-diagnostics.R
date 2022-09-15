library(TwoSampleMR)
library(data.table)
library(MendelianRandomization) 
library(writexl)
library(readxl)
library(simex)
library(MRPRESSO)
library(mr.raps)
library(stringr)
library(cowplot)



exposure_trait = 'CUD'
outcome_trait = '2016EA'

dat = data.frame(fread(paste0('./GWAS/IVS/', exposure_trait, '_to_', outcome_trait, '_instruments.csv'))) 

## Scatter plot:
res <- mr(dat, method_list=c('mr_simple_median', 'mr_weighted_median', 'mr_weighted_mode', 'mr_ivw'))
s <- mr_scatter_plot(res, dat)
# #Funnel plot:
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw"))
fu <- mr_funnel_plot(res_single)
# Forest plot
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw"))
fo <- mr_forest_plot(res_single)
## Leave-one-out plot:
res_loo <- mr_leaveoneout(dat)
l <- mr_leaveoneout_plot(res_loo)


# for 2016EA to LCU
new_plot = plot_grid(s[[1]],fu[[1]],fo[[1]],l[[1]],
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot(paste0("./GWAS/Scatter_plot_",
                 exposure_trait, '_to_', outcome_trait, '.jpg'),
          new_plot,
          ncol = 2, 
          nrow = 2, 
          base_height = 6,
          # base_width = 6, 
          # dpi = 800)

# for 2016EA to CUD
new_plot = plot_grid(s[[1]],fu[[1]],fo[[1]],l[[1]],
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot(paste0("./GWAS/Scatter_plot_",
                 exposure_trait, '_to_', outcome_trait, '.jpg'),
          new_plot,
          ncol = 2, 
          nrow = 2, 
          base_height = 6,
          base_width = 6, 
          dpi = 800)

# for CUD to 2016EA
new_plot = plot_grid(s[[1]],fu[[1]],fo[[1]],l[[1]],
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot(paste0("./GWAS/Scatter_plot_", 
                 exposure_trait, '_to_', outcome_trait, '.pdf'),
          new_plot,
          ncol = 2, 
          nrow = 2, 
          base_height = 3.71,
          base_width = 3.8, 
          dpi = 800)

# for LCU to 2016EA
new_plot = plot_grid(s[[1]],fu[[1]],fo[[1]],l[[1]],
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot(paste0("./GWAS/Scatter_plot_",
                 exposure_trait, '_to_', outcome_trait, '.pdf'),
          new_plot,
          ncol = 2, 
          nrow = 2, 
          base_height = 3.71,
          base_width = 3.8,  
          dpi = 800)


























