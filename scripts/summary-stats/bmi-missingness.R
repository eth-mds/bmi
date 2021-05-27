library(ricu)
library(assertthat)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "utils")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

bmi_proportion <- function(src) {
  
  total <- config("cohort")[[src]][["all"]]
  
  bmi <- config("cohort")[[src]][["bmi"]]
  
  paste0(
    round(
      100 * length(bmi) / length(total),
      1
    ),
    "%"
  )
  
}

hypo_summary <- function(source, patient_ids = config("cohort")[[source]][["bmi"]]) {
  
  ins_hypo <- hypo(source, patient_ids, upto = hours(Inf))
  
  return(
    sprintf("%d (%.1f%%)", 
      nrow(ins_hypo), 100 * nrow(ins_hypo) / length(patient_ids)
    )
  )
}

num_measures <- function(cnc, src, patient_ids) {
  
  if (length(cnc) > 1) return(sum(sapply(cnc, num_measures, src = src, patient_ids = patient_ids)))
  
  nrow(load_concepts(cnc, src, patient_ids = patient_ids, verbose = F))
  
}

db <- c("mimic", "eicu", "hirid", "aumc")

# BMI missingness

print("BMI missingness")

paste0(
  sapply(db, bmi_proportion),
  collapse = ", "
)

# cohort size

print("Cohort size")

paste0(
  sapply(db, function(src) length(config("cohort")[[src]][["bmi"]])),
  collapse = ", "
)

# glucose measurements

print("Glucose:")

paste0(
  sapply(db, function(src) num_measures("glu", src, config("cohort")[[src]][["bmi"]])),
  collapse = ", "
)

# lactate measurements

print("Lactate:")

paste0(
  sapply(db, function(src) num_measures("lact", src, config("cohort")[[src]][["insulin"]])),
  collapse = ", "
)

# liver measurements

print("Liver enzymes:")

paste0(
  sapply(db, function(src) num_measures(c("ast", "alt", "bili"), src, config("cohort")[[src]][["insulin"]])),
  collapse = ", "
)

# MAP measurements

print("MAP:")

paste0(
  sapply(db, function(src) num_measures("map", src, config("cohort")[[src]][["insulin"]])),
  collapse = ", "
)

# insulin measurements

print("Hours of insulin:")

paste0(
  sapply(db, function(src) num_measures("ins", src, config("cohort")[[src]][["insulin"]])),
  collapse = ", "
)

# Hypoglycemia information

print("Hypoglycemia information:")

paste0(
  sapply(db, hypo_summary),
  collapse = ", "
)


