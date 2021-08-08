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

### BMI dose-response

src <- c("aumc", "hirid", "mimic", "eicu")

figures <- list(
  tw_avg_glucose = list(
    title = "Time-weighted average glucose by dataset",
    y_label = "Time-weighted average glucose (mg/dL)",
    coh = "bmi",
    pos.x = 0.2,
    pos.y = 0.8
  ),
  hypo = list(
    title = "Hypoglycemia prevalence by dataset",
    y_label = "Hypoglycemia (proportion)",
    coh = "bmi",
    pos.x = 0.5,
    pos.y = 0.8
  ),
  hypo = list(
    title = "Hypoglycemia in highly monitored group by dataset",
    y_label = "Hypoglycemia (proportion)",
    coh = "bmi",
    subset_fn = high_freq,
    pos.x = 0.5,
    pos.y = 0.8
  ),
  death = list(
    title = "Mortality by dataset",
    y_label = "Mortality (proportion)",
    coh = "bmi",
    pos.x = 0.7,
    pos.y = 0.8
  ),
  death = list(
    title = "Mortality in hypoglycemic group by dataset",
    y_label = "Mortality (proportion)",
    coh = "bmi",
    subset_fn = hypo_only,
    pos.x = 0.7,
    pos.y = 0.8
  ),
  mort_cond_age = list(
    title = "Mortality by age quartile",
    y_label = "Mortality (proportion)",
    coh = "bmi",
    z_cond = "age",
    z_cond_name = "Age Quartile",
    z_cond_labels = paste0("Q", 1:4),
    pos.x = 0.6,
    pos.y = 0.8
  ),
  glu_freq = list(
    title = "Frequency of glucose measurements by dataset",
    y_label = "Mean time between gluc. measurements (h)",
    coh = "bmi",
    pos.x = 0.4,
    pos.y = 0.5
  ),
  fhm = list(
    title = "Value of first hypoglycemic measurement by dataset",
    y_label = "Blood glucose (mg/dL)",
    coh = "bmi",
    pos.x = 0.4,
    pos.y = 0.8
  )#,
  # hypo_cond_diab = list(
  #   title = "Hypoglycemia and diabetes",
  #   y_label = "Hypoglycemia (proportion)",
  #   coh = "bmi",
  #   z_cond = "DM",
  #   z_cond_name = "Diabetes",
  #   z_cond_labels = c("No", "Yes"),
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # death_cond_diab = list(
  #   title = "Mortality and diabetes",
  #   y_label = "Mortality (proportion)",
  #   coh = "bmi",
  #   z_cond = "DM",
  #   z_cond_name = "Diabetes",
  #   z_cond_labels = c("No", "Yes"),
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # hypo_cond_adm = list(
  #   title = "Hypoglycemia and admission type",
  #   y_label = "Hypoglycemia (proportion)",
  #   coh = "bmi",
  #   z_cond = "adm",
  #   z_cond_name = "Admission Type",
  #   z_cond_labels = c("Medical", "Surgical", "Other"),
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # death_cond_adm = list(
  #   title = "Mortality and admission type",
  #   y_label = "Mortality (proportion)",
  #   coh = "bmi",
  #   z_cond = "adm",
  #   z_cond_name = "Admission Type",
  #   z_cond_labels = c("Medical", "Surgical", "Other"),
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # max_insulin = list(
  #   title = "Insulin rate by dataset",
  #   y_label = "Maximal insulin rate (U/h)",
  #   coh = "insulin",
  #   pos.x = 0.2,
  #   pos.y = 0.8
  # ),
  # max_insulin_wnorm = list(
  #   title = "Weight-normalized insulin by dataset",
  #   y_label = "Weight-norm. max. insulin rate (U/(h*kg))",
  #   coh = "insulin",
  #   pos.x = 0.3,
  #   pos.y = 0.8
  # ),
  # dur_TPN = list(
  #   title = "Total parenteral nutrition duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # dur_enteral = list(
  #   title = "Enteral nutrition duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # dur_cortico = list(
  #   title = "Corticosteroids duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # ),
  # tw_avg_dextrose = list(
  #   title = "Average dextrose infusion rate",
  #   y_label = "Infusion rate (ml/hour)",
  #   coh = "insulin",
  #   pos.x = 0.6,
  #   pos.y = 0.8
  # )
)

# conditional, multi-source plots
age_binning <- function(x)
  .bincode(x, c(-Inf, quantile(x, c(0.25, 0.5, 0.75), na.rm = FALSE), Inf))

diab_src <- c("mimic", "eicu")
binary_binning <- function(x) .bincode(x, c(-Inf, 0.5, Inf))

figures[["hypo_cond_diab"]][["dat"]] <- CI_dat(diab_src, y = "hypo", z = "DM",
                                               z_binning = binary_binning)
figures[["death_cond_diab"]][["dat"]] <- CI_dat(diab_src, y = "death", z = "DM",
                                                z_binning = binary_binning)

adm_src <- c("aumc", "mimic", "eicu")
adm_binning <- function(x)
  ifelse(x == "med", 1, ifelse(x == "surg", 2, NA))
figures[["hypo_cond_adm"]][["dat"]] <- CI_dat(adm_src, y = "hypo", z = "adm",
                                              z_binning = adm_binning)
figures[["death_cond_adm"]][["dat"]] <- CI_dat(adm_src, y = "death", z = "adm",
                                               z_binning = adm_binning)

figures[["mort_cond_age"]][["dat"]] <- CI_dat(src, y = "death", z = "age",
                                              z_binning = age_binning)

# unconditional, single-source plots
for (i in 1:length(figures)) {
  
  if (is.null(figures[[i]][["dat"]])) {
    
    cat("Handling target", names(figures)[i], "\n\n")
    dat <- lapply(src, CI_dat, y = names(figures)[i], 
                  subset_fn = figures[[i]][["subset_fn"]])
    figures[[i]][["dat"]] <- dat
    
  } else dat <- figures[[i]][["dat"]]
  
  cat("-----------------------\n")
  
  plot <- CI_plot(dat, 
                  title = figures[[i]][["title"]], 
                  y_label = figures[[i]][["y_label"]],
                  z_cond = figures[[i]][["z_cond"]],
                  z_cond_name = figures[[i]][["z_cond_name"]],
                  z_cond_labels = figures[[i]][["z_cond_labels"]],
                  pos.x = figures[[i]][["pos.x"]],
                  pos.y = figures[[i]][["pos.y"]])
  
  
  figures[[i]][["plot"]] <- plot
  
}

# Figures
fig1 <- plot_grid(
  figures[["tw_avg_glucose"]][["plot"]],
  figures[["death"]][["plot"]],
  figures[["mort_cond_age"]][["plot"]],
  figures[["hypo"]][["plot"]],
  labels = c("a)", "b)", "c)", "d)"), 
  ncol = 2L
)

ggsave(file.path(root, "figures", "Figure3.tiff"),
       plot = fig1, width = 18, height = 12)

fig2 <- plot_grid(
  figures[["death_cond_adm"]][["plot"]],
  figures[["death_cond_diab"]][["plot"]],
  figures[["hypo_cond_adm"]][["plot"]],
  figures[["hypo_cond_diab"]][["plot"]],
  labels = c("a)", "b)", "c)", "d)"), 
  ncol = 2L
)

ggsave(file.path(root, "figures", "Figure4.tiff"),
       plot = fig2, width = 18, height = 12)

leg3 <- get_legend(
  figures[["glu_freq"]][["plot"]] + theme(legend.position = "bottom")
)
fig3 <- plot_grid(
  figures[["glu_freq"]][["plot"]] + theme(legend.position = "none"),
  figures[[3L]][["plot"]] + theme(legend.position = "none"),
  figures[["fhm"]][["plot"]] + theme(legend.position = "none"),
  figures[[5L]][["plot"]] + theme(legend.position = "none"),
  labels = c("a)", "b)", "c)", "d)"), 
  ncol = 2L
)
fig3 <- plot_grid(fig3, leg3, ncol = 1L, rel_heights = c(1, 0.075))

ggsave(file.path(root, "figures", "Figure5.tiff"),
       plot = fig3, width = 18, height = 12)

leg4 <- get_legend(
  figures[["max_insulin"]][["plot"]] + theme(legend.position = "bottom")
)
fig4 <- plot_grid(
  figures[["max_insulin"]][["plot"]] + theme(legend.position = "none"),
  figures[["max_insulin_wnorm"]][["plot"]] + theme(legend.position = "none"),
  figures[["dur_TPN"]][["plot"]] + theme(legend.position = "none"),
  figures[["dur_enteral"]][["plot"]] + theme(legend.position = "none"),
  figures[["dur_cortico"]][["plot"]] + theme(legend.position = "none"),
  figures[["tw_avg_dextrose"]][["plot"]] + theme(legend.position = "none"),
  labels = c("a)", "b)", "c)", "d)", "e)", "f)"), 
  ncol = 2L
)
fig4 <- cowplot::plot_grid(fig4, leg4, ncol = 1L, rel_heights = c(1, 0.05))

ggsave(file.path(root, "figures", "Figure6.tiff"),
       plot = fig4, width = 18, height = 18)
