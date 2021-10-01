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
  gv_cv = list(
    title = "Glucose coefficient of variation by dataset",
    y_label = "Blood glucose coefficient of variation (%)",
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
  glu_freq = list(
    title = "Frequency of glucose measurements by dataset",
    y_label = "Mean time between gluc. measurements (h)",
    coh = "bmi",
    pos.x = 0.4,
    pos.y = 0.5
  ),
  hypo_dur = list(
    title = "Duration of hypoglycemia by dataset",
    y_label = "Proportion of glucose values below 70 mg/dL",
    coh = "bmi",
    subset_fn = hypo_only,
    pos.x = 0.4,
    pos.y = 0.8
  ),
  hypo_cnt = list(
    title = "Number of hypoglycemic episodes by dataset",
    y_label = "Number of hypoglycemic episodes",
    coh = "bmi",
    subset_fn = hypo_only,
    pos.x = 0.4,
    pos.y = 0.8
  ),
  hypo_sev = list(
    title = "Lowest glucose level by dataset",
    y_label = "Blood glucose (mg/dL)",
    coh = "bmi",
    subset_fn = hypo_only,
    pos.x = 0.4,
    pos.y = 0.8
  )
  # max_insulin = list(
  #   title = "Insulin rate by dataset",
  #   y_label = "Maximal insulin rate (U/h)",
  #   coh = "insulin",
  #   pos.x = 0.2,
  #   pos.y = 0.8, add_prop = TRUE
  # ),
  # max_insulin_wnorm = list(
  #   title = "Weight-normalized insulin by dataset",
  #   y_label = "Weight-norm. max. insulin rate (U/(h*kg))",
  #   coh = "insulin",
  #   pos.x = 0.3,
  #   pos.y = 0.8, add_prop = TRUE
  # ),
  # dur_TPN = list(
  #   title = "Total parenteral nutrition duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.6,
  #   pos.y = 0.8, add_prop = TRUE
  # ),
  # dur_enteral = list(
  #   title = "Enteral nutrition duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.5,
  #   pos.y = 0.6, add_prop = TRUE
  # ),
  # dur_cortico = list(
  #   title = "Corticosteroids duration",
  #   y_label = "Mean proportion of time treated",
  #   coh = "insulin",
  #   pos.x = 0.5,
  #   pos.y = 0.8, add_prop = TRUE
  # ),
  # tw_avg_dextrose = list(
  #   title = "Average dextrose infusion rate",
  #   y_label = "Infusion rate (ml/hour)",
  #   coh = "insulin",
  #   pos.x = 0.35,
  #   pos.y = 0.8, add_prop = TRUE
  # )
)

# unconditional, single-source plots

for (i in c(9, 10)) {
#for (i in 1:length(figures)) {
  
  cat("-----------------------\n")
  
  cat("Handling target", names(figures)[i], "\n\n")
  dat <- lapply(src, CI_dat, y = names(figures)[i],
                subset_fn = figures[[i]][["subset_fn"]],
                add_prop = figures[[i]][["add_prop"]])
  figures[[i]][["dat"]] <- dat
  
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

# leg1 <- get_legend(
#   figures[["tw_avg_glucose"]][["plot"]] + theme(legend.position = "bottom")
# )
# 
# fig1 <- plot_grid(
#   figures[["tw_avg_glucose"]][["plot"]] +theme(legend.position = "none"),
#   figures[["death"]][["plot"]] + theme(legend.position = "none"),
#   figures[["hypo"]][["plot"]] + theme(legend.position = "none"),
#   figures[["gv_cv"]][["plot"]] + theme(legend.position = "none"),
#   labels = c("a)", "b)", "c)", "d)"), 
#   ncol = 2L
# )
# fig1 <- plot_grid(fig1, leg1, ncol = 1L, rel_heights = c(1, 0.05))
# 
# ggsave(file.path(root, "figures", "Figure2.0.tiff"), plot = fig1, 
#        width = 18, height = 12, type = "cairo", compression = "lzw")
# 
leg3 <- get_legend(
  figures[["glu_freq"]][["plot"]] + theme(legend.position = "bottom")
)

fig3 <- plot_grid(
  figures[["glu_freq"]][["plot"]] + theme(legend.position = "none"),
  figures[[4L]][["plot"]] + theme(legend.position = "none"),
  figures[[6L]][["plot"]] + theme(legend.position = "none"),
  figures[["hypo_sev"]][["plot"]] + theme(legend.position = "none"),
  figures[["hypo_dur"]][["plot"]] + theme(legend.position = "none"),
  figures[["hypo_cnt"]][["plot"]] + theme(legend.position = "none"),
  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
  ncol = 2L
)
fig3 <- plot_grid(fig3, leg3, ncol = 1L, rel_heights = c(1, 0.05))

ggsave(file.path(root, "figures", "Figure5.0.tiff"), plot = fig3,
       width = 18, height = 18, type = "cairo", compression = "lzw")

# fig4 <- plot_grid(
#   figures[["max_insulin"]][["plot"]],
#   figures[["max_insulin_wnorm"]][["plot"]],
#   figures[["dur_TPN"]][["plot"]],
#   figures[["dur_enteral"]][["plot"]],
#   figures[["dur_cortico"]][["plot"]],
#   figures[["tw_avg_dextrose"]][["plot"]],
#   labels = c("a)", "b)", "c)", "d)", "e)", "f)"), 
#   ncol = 2L
# )
# 
# ggsave(file.path(root, "figures", "Figure6.0.tiff"), plot = fig4, 
#        width = 18, height = 18, type = "cairo", compression = "lzw")
