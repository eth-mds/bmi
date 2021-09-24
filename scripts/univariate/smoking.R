library(ricu)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(boot)
library(data.table)
library(icd)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

smk <- CI_dat("mimic", y = "SMK")
pfi <- lapply(c("mimic", "eicu", "hirid", "aumc"), CI_dat, y = "min_pafi")
plot_grid(
  CI_plot(list(smk), y_label = "Proportion smokers", 
          title = "Smoking MIMIC-III"),
  CI_plot(pfi, y_label = "Worst PaO2/FiO2 value", 
          title = "PaO2/FiO2 at 24 hours"),
  ncol = 2L
)

ggsave(file.path(root, "figures", "smoking.tiff"), height = 6, width = 12,
       compression = "lzw", type = "cairo")
