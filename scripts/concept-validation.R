library(assertthat)
library(ricu)
library(ggplot2)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

conc_val <- function(cnc, src) {
  
  coh <- config("cohort")[[src]][["insulin"]]
  
  tbl <- load_concepts(cnc, src, patient_ids = coh)
  prop <- length(unique(id_col(tbl))) / length(coh)
  
  # need to fix the duration of treatment -> using expand()!
  
  tbl <- expand(tbl, aggregate = TRUE)
  tbl[, list(total = sum(get(cnc))), by = c(id_vars(tbl))]
  tbl <- merge(tbl, load_concepts("los_icu", src), all.x = TRUE)
   
  tbl[, total := total / los_icu]
  c(median(tbl$total), prop)
  
}

res <- NULL
src <- c("mimic", "eicu", "hirid", "aumc")
for (c in c("TPN", "enteral", "cortico")) {
  
  for (dsrc in src) {
    
    res <- rbind(
      res, 
      data.frame(conc_val(c, dsrc), concept = c, source = dsrc)
    )   
    
  }
  
}

ggplot(res, aes(x = concept, y = prop, fill = source)) +
  geom_col(position = "dodge") + theme_bw() +
  xlab("Concept") + ylab("% treated") + ggtitle("Proportion of treated patients")

# example...

# tbl <- load_concepts("cortico", "aumc")
# tbl[, 
#     c(dur_var(tbl)) := as.difftime(ceiling(as.numeric(get(dur_var(tbl)) / 60)), 
#                                      units = "hours")]
# tbl[, end_var := get(index_var(tbl)) + get(dur_var(tbl))]
# 
# tbl <- as_ts_tbl(tbl)
# 
# expand(tbl, end_var = "end_var", aggregate = "min", keep_vars = "cortico")
# 
