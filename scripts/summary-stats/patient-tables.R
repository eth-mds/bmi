library(ricu)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(plyr)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "utils")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- c("aumc", "hirid", "mimic", "eicu")
cohorts <- lapply(src, function(x) config("cohort")[[x]][["bmi"]])

vars <- list(
  age = list(
    concept = "age",
    callback = med_iqr
  ),
  bmi = list(
    concept = "bmi_bins",
    callback = tab_design
  ),
  admission = list(
    concept = "adm",
    callback = tab_design
  ),
  death = list(
    concept = "death",
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

res <- Reduce(
  function(x, y) merge(x, y, by = c("Variable", "Reported"), sort = F),
  Map(pts_source_sum, src, cohorts)
)

my_doc <- read_docx()

my_doc <- my_doc %>%
  body_add_table(res, style = "table_template")

print(my_doc, target = file.path(root, "tables", "Table1.docx"))
