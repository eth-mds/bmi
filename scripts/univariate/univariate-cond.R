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

cndl <- list(
  diabetes = list(
    src = c("mimic", "eicu"),
    z = "DM",
    z_binning = function(x) .bincode(x, c(-Inf, 0.5, Inf)),
    titles = "by diabetes status",
    z_cond_name = "Diabetes",
    z_cond_labels = c("No", "Yes"),
    figname = "Figure3.0"
  ),
  hba1c = list(
    src = "mimic",
    z = "hba1c",
    z_binning = function(x) .bincode(x, c(-Inf, 6.1, 6.5, 7, Inf)),
    titles = "by HbA1c level",
    z_cond_name = "HbA1c",
    z_cond_labels = c("<6.1%", "6.1-6.5%", "6.6-7%", ">7%"),
    figname = "Figure3.5"
  )#,
  # admission = list(
  #   src = c("aumc", "mimic", "eicu"),
  #   z = "adm",
  #   z_binning = function(x) ifelse(x == "med", 1, ifelse(x == "surg", 2, NA)),
  #   titles = "by admission type",
  #   z_cond_name = "Admission",
  #   z_cond_labels = c("Medical", "Surgical", "Other"),
  #   figname = "Figure2.5"
  # )
)

out <- c("tw_avg_glu", "death", "hypo", "gv_cv")
ttl <- c("Time-weighted average glucose",
         "Mortality", "Hypoglycemia", "Glucose variability")
ylabs <- c("Time-weighted glucose average (mg/dL)",
           "Mortality (proportion)", "Hypoglycemia (proportion)",
           "Blood glucose coefficient of variation (%)")

for (i in seq_along(cndl)) {
  
  plt <- list()
  for (j in seq_along(out)) {
    
    res <- CI_dat(cndl[[i]][["src"]], y = out[j], z = cndl[[i]][["z"]],
                  z_binning = cndl[[i]][["z_binning"]])
    plt[[j]] <- CI_plot(
      res, title = paste(ttl[j], cndl[[i]][["titles"]]),
      y_label = ylabs[j], z_cond = cndl[[i]][["z"]],
      z_cond_name = cndl[[i]][["z_cond_name"]], 
      z_cond_labels = cndl[[i]][["z_cond_labels"]]
    )
    leg <- get_legend(plt[[j]] + theme(legend.position = "bottom"))
    plt[[j]] <- plt[[j]] + theme(legend.position = "none")
    cat("\n\n", i, "and ", out[j], "\n\n")
    
  }
  fig <- plot_grid(plotlist = plt, ncol = 2L, labels = c("a)", "b)", "c)", "d)"))
  fig <- plot_grid(fig, leg, ncol = 1L, rel_heights = c(1, 0.05))
  ggsave(file.path(root, "figures", paste0(cndl[[i]][["figname"]], ".tiff")),
         plot = fig, width = 18, height = 12, compression = "lzw", type = "cairo")
  
}
