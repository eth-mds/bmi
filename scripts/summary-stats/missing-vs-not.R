library(ricu)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(plyr)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

vars <- list(
  age = list(
    concept = "age",
    callback = med_iqr
  ),
  admission = list(
    concept = "adm",
    callback = tab_design
  ),
  death = list(
    concept = "death",
    callback = percent_fun
  ),
  is_hypo = list(
    concept = "is_hypo",
    callback = percent_fun
  ),
  los_icu = list(
    concept = "los_icu",
    callback = med_iqr
  ),
  los_hosp = list(
    concept = "los_hosp",
    callback = med_iqr
  ),
  gender = list(
    concept = "sex",
    callback = tab_design
  ),
  sofa = list(
    concept = "sofa",
    callback = multi_med_iqr
  )
)

all <- lapply(config("cohort"), `[[`, "all")
bmi <- lapply(config("cohort"), `[[`, "bmi")
miss <- Map(function(x, y) setdiff(x, y), all, bmi)

src <- c("mimic_demo", "eicu_demo")
names(miss) <- c(src, "hirid", "aumc")

pts_source_sum(src, patient_ids = miss[1:2])
