library(ricu)
library(assertthat)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

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

coh_size <- function(src, coh) {
  Reduce(sum, lapply(config("cohort")[src], function(x) length(x[[coh]])))
}

hypo_summary <- function(src, coh) {
  
  patient_ids <- config("cohort")[[src]][[coh]]
  
  ins_hypo <- hypo(source, patient_ids, upto = hours(Inf))
  
  return(
    sprintf("%d (%.1f%%)", 
      nrow(ins_hypo), 100 * nrow(ins_hypo) / length(patient_ids)
    )
  )
}

num_measures <- function(cnc, src, coh) {
  
  patient_ids <- lapply(config("cohort")[src], `[[`, coh)
  tbl <- load_concepts(cnc, src, patient_ids = patient_ids, verbose = F)
  cat("Concept", cnc, "with a total of", nrow(tbl), "measures\n")
  if (length(src) > 1L) 
    cat(paste(tbl[, .N, by = "source"][["N"]], collapse = ", "), "per dataset\n")
}

db <- c("aumc", "hirid", "mimic", "eicu")

# BMI missingness

print("BMI missingness")
cat("overall:", 100 * coh_size(db, "bmi") / coh_size(db, "all"))
paste0(
  sapply(db, bmi_proportion),
  collapse = ", "
)

# cohort size
print("Cohort size")
bmi <- sapply(db, coh_size, coh = "bmi")
sum(bmi)
paste0(bmi, collapse = ", ")

# glucose measurements
num_measures("glu", db, "bmi")

# multivariate analysis
mcnc <- c("lact", "map", "norepi_equiv", "ins_ifx", "dex_amount", "TPN", 
          "enteral", "cortico")
for (cnc in mcnc) num_measures(cnc, db, "insulin")

# Hypoglycemia information
print("Hypoglycemia information:")

paste0(
  sapply(db, hypo_summary),
  collapse = ", "
)


