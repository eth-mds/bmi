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

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

### BMI dose-response

src <- c("mimic", "eicu", "hirid", "aumc")

figures <- list(
  tw_avg_glucose = list(
    title = "Time-weighted average glucose by dataset",
    y_label = "Time-weighted average glucose (mg/dL)",
    coh = "bmi",
    pos.x = 0.2,
    pos.y = 0.8
  ),
  death = list(
    title = "Mortality by dataset",
    y_label = "Mortality (proportion)",
    coh = "bmi",
    pos.x = 0.7,
    pos.y = 0.8
  ),
  hypo = list(
    title = "Hypoglycemia prevalence by dataset",
    y_label = "Hypoglycemic propensity (proportion)",
    coh = "bmi",
    pos.x = 0.5,
    pos.y = 0.8
  ),
  max_insulin = list(
    title = "Insulin rate by dataset",
    y_label = "Maximal insulin rate (U/h)",
    coh = "insulin",
    pos.x = 0.2,
    pos.y = 0.8
  ),
  max_insulin_wnorm = list(
    title = "Weight-normalized insulin by dataset",
    y_label = "Weight-norm. max. insulin rate (U/(h*kg))",
    coh = "insulin",
    pos.x = 0.3,
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
  ),
  hypo = list(
    title = "Hypoglycemia in highly monitored group by dataset",
    y_label = "Hypoglycemic propensity (proportion)",
    coh = "bmi",
    subset_fn = high_freq,
    pos.x = 0.5,
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
  twavg_cond_death = list(
    title = "Time-weighted glucose by in-hospital survival",
    y_label = "Time-weighted glucose (mg/dL)",
    coh = "bmi",
    z_cond = "death",
    z_cond_name = "In-hospital death",
    z_cond_labels = c("FALSE", "TRUE"),
    pos.x = 0.2,
    pos.y = 0.8
  )
)

figures <- list(
  hypo_cond_diab = list(
    title = "Hypoglycemia for Diabetic vs. Non-diabetic",
    y_label = "Hypoglycemia (proportion)",
    coh = "bmi",
    z_cond = "DM",
    z_cond_name = "Diabetes",
    z_cond_labels = c("No", "Yes"),
    pos.x = 0.6,
    pos.y = 0.8
  ),
  death_cond_diab = list(
    title = "Death for Diabetic vs. Non-diabetic",
    y_label = "Death (proportion)",
    coh = "bmi",
    z_cond = "DM",
    z_cond_name = "Diabetes",
    z_cond_labels = c("No", "Yes"),
    pos.x = 0.6,
    pos.y = 0.8
  ),
  hypo_cond_adm = list(
    title = "Hypoglycemia and admission type",
    y_label = "Hypoglycemia (proportion)",
    coh = "bmi",
    z_cond = "adm",
    z_cond_name = "Admission Type",
    z_cond_labels = c("Medical", "Surgical", "Other"),
    pos.x = 0.6,
    pos.y = 0.8
  ),
  death_cond_adm = list(
    title = "Death and admission type",
    y_label = "Death (proportion)",
    coh = "bmi",
    z_cond = "adm",
    z_cond_name = "Admission Type",
    z_cond_labels = c("Medical", "Surgical", "Other"),
    pos.x = 0.6,
    pos.y = 0.8
  )
)

# conditional, multi-source plots
age_binning <- function(x) 
  .bincode(x, c(-Inf, quantile(x, c(0.25, 0.5, 0.75), na.rm = FALSE), Inf))

bin_binning <- function(x) .bincode(x, c(-Inf, 0.5, Inf))

adm_binning <- function(x) 
  ifelse(x == "med", 1, ifelse(x == "surg", 2, NA))

figures[["hypo_cond_diab"]][["dat"]] <- CI_dat(c("mimic"), y = "hypo", 
                                               z = "DM", z_binning = bin_binning)
figures[["death_cond_diab"]][["dat"]] <- CI_dat(c("mimic"), y = "death", 
                                               z = "DM", z_binning = bin_binning)

figures[["hypo_cond_adm"]][["dat"]] <- CI_dat(c("mimic"), y = "hypo", 
                                               z = "adm", z_binning = adm_binning)
figures[["death_cond_adm"]][["dat"]] <- CI_dat(c("mimic"), y = "death", 
                                                z = "adm", z_binning = adm_binning)

figures[["mort_cond_age"]][["dat"]] <- CI_dat(src, y = "death", z = "age",
                                              z_binning = age_binning)

figures[["twavg_cond_death"]][["dat"]] <- CI_dat(src, 
                                                 y = "tw_avg_glucose", 
                                                 z = "death",
                                                 z_binning = function(x) x)

# unconditional, single-source plots
for (i in 1:length(figures)) {
  
  if (is.null(figures[[i]][["dat"]])) {
    
    dat <- lapply(src, CI_dat, y = names(figures)[i])
    figures[[i]][["dat"]] <- dat
    
  } else dat <- figures[[i]][["dat"]]
  
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
Figure_2 <- cowplot::plot_grid(
  figures[["tw_avg_glucose"]][["plot"]],
  figures[["death"]][["plot"]],
  figures[["mort_cond_age"]][["plot"]],
  figures[["hypo"]][["plot"]], 
  labels = c("a)", "b)", "c)", "d)"), 
  ncol = 2L
)

# ggsave(file.path(root, "4files-BMI", "figures", "Figure2.tiff"), 
#        plot = Figure_2, width = 18, height = 12)

eFigure_1 <- cowplot::plot_grid(
  figures[[8L]][["plot"]],
  figures[["twavg_cond_death"]][["plot"]],
  figures[["max_insulin"]][["plot"]],
  figures[["max_insulin_wnorm"]][["plot"]], labels = c("a)", "b)", "c)", "d)"), 
  ncol = 2L
)

# ggsave(file.path(root, "4files-BMI", "eSupplement", "eFigure1.tiff"), 
#        plot = eFigure_1, width = 18, height = 12)

eFigure_2 <- cowplot::plot_grid(
  figures[["glu_freq"]][["plot"]],
  figures[[9L]][["plot"]],
  figures[["fhm"]][["plot"]],
  labels = c("a)", "b)", "c)"), 
  ncol = 3L
)

# ggsave(file.path(root, "4files-BMI", "eSupplement", "eFigure2.tiff"), 
#        plot = eFigure_2, width = 24, height = 6)
