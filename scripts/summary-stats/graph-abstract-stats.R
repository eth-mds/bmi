library(ricu)
library(assertthat)
library(boot)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- c("aumc", "hirid", "mimic", "eicu")
bmi <- lapply(config("cohort"), `[[`, "bmi")
y_range <- c("death", "hypo", "tw_avg_glucose")

res <- NULL
for (y in y_range) {
  
  ybin <- CI_dat(src, y = y, x_bins = c(18.5, 25), patient_ids = bmi)
  ybin <- data.table::setnames(ybin, c("meanval", "V1"), c("Group", "Value"))
  if (is.element(y, c("death", "hypo"))) {
    ybin[, Value := 100 * Value]
  }
  ybin <- ybin[, c("Group", "Value"), with = FALSE]
  ybin[, Feature := y]
  
  res <- rbind(res, ybin)
  
}

res
