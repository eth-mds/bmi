library(assertthat)
library(ricu)
library(ggplot2)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

conc_val <- function(cnc, src) {
  
  coh <- config("cohort")[[src]][["insulin"]]
  
  tbl <- load_concepts(cnc, src)
  tbl <- tbl[get(id_var(tbl)) %in% coh]
  prop <- length(unique(id_col(tbl))) / length(coh)
  
  # need to fix the duration of treatment -> using expand()!
  tbl <- expand(tbl, aggregate = TRUE)
  tbl <- tbl[, list(total = .N), by = c(id_vars(tbl))]
  tbl <- merge(tbl, load_concepts("los_icu", src), all.x = TRUE)
   
  tbl[, total := total / los_icu]
  c(median(tbl$total), prop)
  
}

res <- NULL
src <- c("mimic", "eicu", "hirid", "aumc")
for (c in c("TPN", "enteral", "cortico")) {
  
  for (dsrc in src) {
    
    out <- conc_val(c, dsrc)
    prop <- out[2]
    dpd <- out[1]
    res <- rbind(
      res, 
      data.frame(prop = prop, dpd = dpd, concept = c, source = dsrc)
    )   
    
  }
  
}

prop <- ggplot(res, aes(x = concept, y = prop, fill = source)) +
  geom_col(position = "dodge") + theme_bw() +
  xlab("Concept") + ylab("% treated") + ggtitle("Proportion of treated patients")

dur <- ggplot(res, aes(x = concept, y = dpd, fill = source)) +
  geom_col(position = "dodge") + theme_bw() +
  xlab("Concept") + ylab("Treatment duration (% of day)") + 
  ggtitle("Median duration of treatment per day of patient stay")

cowplot::plot_grid(prop, dur, ncol = 2L)