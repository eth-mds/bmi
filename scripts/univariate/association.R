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

src <-  c("mimic", "eicu")
pids <- lapply(config("cohort")[src], `[[`, "bmi")
res <- load_concepts(c("bmi_bins", "is_hypo", "DM", "death", "gv_cv",
                       "tw_avg_glu", "hypo_cnt"), src,
                     patient_ids = pids,
                     verbose = FALSE)
res[is.na(death), "death"] <- FALSE
res[is.na(DM), "DM"] <- 0

assoc_hmap(res, "tw_avg_glu", bins = c(140, 180),
           c("70-140 mg/dL", "140-180 mg/dL", "> 180 mg/dL"),
           y_label = "Time-weighted average glucose",
           title = "Average glucose and mortality")
ggsave(file.path(root, "figures", "Figure4.0.tiff"), 
       width = 15, height = 6, compression = "lzw", type = "cairo")

assoc_hmap(res)
ggsave(file.path(root, "figures", "Figure4.2.tiff"), 
       width = 15, height = 6, compression = "lzw", type = "cairo")

assoc_hmap(res, "hypo_cnt", bins = c(0.5, 1.5), c("0", "1", "> 1"),
           y_label = "Number of hypoglycemic episodes",
           title = "Hypoglycemia and mortality")
ggsave(file.path(root, "figures", "Figure4.3.tiff"),
       width = 15, height = 6, compression = "lzw", type = "cairo")

### Non-hypoglycemic only
pids_hypo <- Map(hypo_only, src, pids, hours(Inf))
pids_nhypo <- Map(setdiff, pids, pids_hypo)
res_nhypo <- load_concepts(c("bmi_bins", "is_hypo", "DM", "death", "gv_cv",
                           "tw_avg_glu", "hypo_cnt"), src, 
                           patient_ids = pids_nhypo, verbose = FALSE)
res_nhypo[is.na(death), "death"] <- FALSE
res_nhypo[is.na(DM), "DM"] <- 0
assoc_hmap(res_nhypo, "tw_avg_glu", bins = c(140, 180),
           c("70-140 mg/dL", "140-180 mg/dL", "> 180 mg/dL"),
           y_label = "Time-weighted average glucose",
           title = "Average glucose and mortality (non-hypoglycemic group)")
ggsave(file.path(root, "figures", "Figure4.1.tiff"), 
       width = 15, height = 6, compression = "lzw", type = "cairo")
