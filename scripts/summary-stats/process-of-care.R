library(ricu)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(plyr)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "utils")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

src <- c("mimic", "eicu", "hirid", "aumc")
cohorts <- lapply(src, function(x) config("cohort")[[undemo(x)]][["insulin"]])

prop_fun <- function(var_name, var, patient_ids) {
  list(var_name, "%", round(100*sum(var, na.rm = T) / length(patient_ids)))
}

med_iqr_fun <- function(var_name, var, patient_ids) {
  mlp <- paste0(
    round(quantile(var, 0.5, na.rm = T), 2), " (",
    round(quantile(var, 0.25, na.rm = T), 2), "-",
    round(quantile(var, 0.75, na.rm = T), 2), ")"
  )
  list(var_name, "Median (IQR)", mlp)
}

vars <- list(
  hypo = list(
    concept = "glu",
    direction = "decreasing",
    threshold = 3.9 * 18.016,
    name = "Proportion with hypoglycemia",
    callback = prop_fun
  )#,
  # insulin = list(
  #   concept = "ins",
  #   direction = "increasing",
  #   threshold = 0,
  #   name = "Proportion on insulin",
  #   callback = prop_fun
  # ),
  # vasopressors = list(
  #   concept = "norepi_equiv",
  #   direction = "increasing",
  #   threshold = 0,
  #   name = "Proportion on vasopressors",
  #   callback = prop_fun
  # ),
  # lactate = list(
  #   concept = "lact",
  #   direction = "increasing",
  #   name = "Maximal lactate (mmol/L)",
  #   callback = med_iqr_fun
  # ),
  # bili = list(
  #   concept = "bili",
  #   direction = "increasing",
  #   name = "Maximal bilirubin (mg/dL)",
  #   callback = med_iqr_fun
  # ),
  # ast = list(
  #   concept = "ast",
  #   direction = "increasing",
  #   name = "Maximal AST (IU/L)",
  #   callback = med_iqr_fun
  # ),
  # alt = list(
  #   concept = "alt",
  #   direction = "increasing",
  #   name = "Maximal ALT (IU/L)",
  #   callback = med_iqr_fun
  # ),
  # meanbp = list(
  #   concept = "map",
  #   direction = "decreasing",
  #   name = "Minimal MAP (mmHg)",
  #   callback = med_iqr_fun
  # )
)

gen_cov <- function(src, concept, direction, threshold, patient_ids) {

  if (concept == "tw_avg_glucose") {
      return(rnorm(100))
  } else if (concept == "gluc_meas_pd") {
      return(rnorm(100))
  } else {

    cov <- w_value(src, concept, direction, upto = hours(72L), 
      patient_ids = patient_ids, imp_val = NA_real_)[["w_val"]]

    if(!is.null(threshold)) {

      if (direction != "increasing") return(cov < threshold)
      return(cov > threshold)

    } else return(cov)

  }

}

poc_source_sum <- function(source, patient_ids) {

  tbl_list <- lapply(
    vars,
    function(x) x[["callback"]]( x[["name"]],
        gen_cov(source, x[["concept"]], x[["direction"]], x[["threshold"]], patient_ids), patient_ids
      )
  )

  poc_tbl <- Reduce(rbind,
    lapply(
      tbl_list,
      function(x) data.frame(Reduce(cbind, x))
    )
  )

  names(poc_tbl) <- c("Variable (first 72 hours)", "Reported", srcwrap(source))

  poc_tbl

}

res <- Reduce(
  function(x, y) merge(x, y, by = c("Variable (first 72 hours)", "Reported"), sort = F),
  Map(poc_source_sum, src, cohorts)
)

my_doc <- read_docx()

my_doc <- my_doc %>%
  body_add_table(res, style = "table_template")

print(my_doc, target = file.path(root, "tables", "Table2.docx"))
