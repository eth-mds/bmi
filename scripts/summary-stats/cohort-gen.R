library(ricu)
library(ggplot2)
library(assertthat)

source(file.path(rprojroot::find_root(".gitignore"), "utils", "utils-config.R"))

# filter for age on MIMIC-III
mimic <- load_concepts(c("age", "bmi"), "mimic", verbose = F)

# filter for age on MIMIC-III
eicu <- load_concepts(c("age", "bmi"), "eicu", verbose = F)

# filter for age on HiRID
hirid <- load_concepts(c("age", "bmi"), "hirid", verbose = F)

# filter for age on AUMC
aumc <- load_concepts(c("age", "bmi"), "aumc", verbose = F)

# sort out eICU
filter0 <- list(
  lact = -1,
  glu = -1,
  ins = -1,
  bili = -1,
  ast = -1
)

filter1 <- list(
  lact = 0.4,
  glu = 0,
  ins = 0.1,
  bili = 0,
  ast = 0
)

pop_per_hosp <- function(filter) {
  concepts <- names(filter)
  source <- "eicu"
  tbl <- load_concepts(concepts, source, verbose = FALSE)
  res <- tbl[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols = concepts, 
             by = eval(id_vars(tbl))]
  hospitals <- load_id("patient", "eicu")
  hospitals <- hospitals[, c("patientunitstayid", "hospitalid"), with = FALSE]
  hospitals[, hsize := .N, by = "hospitalid"]
  res <- merge(res, hospitals, by = id_vars(tbl))
  res <- res[, c(lapply(.SD, function(x) sum(x > 0)), hsize = hsize[1]), 
             .SDcols = concepts, by = "hospitalid"]
  res <- res[, c(lapply(.SD, function(x) x / hsize[1]), hsize = hsize[1]), 
             .SDcols = concepts, by = "hospitalid"]

  p <- ggplot(data=res, aes(x=factor(hospitalid), y=lact)) +
    geom_bar(stat="identity") + xlab("Measures per patient") + ylab("Hospital") + 
    theme_minimal() + ggtitle("barplot")

  id <- res[, lapply(concepts, function(x) get(x) > filter[[x]])]
  id <- apply(id, 1L, all)

  hid <- res[id][["hospitalid"]]
  print(paste("Hospitals included:", length(hid)))
  print(paste("Number out:", nrow(hospitals[!(hospitalid %in% hid)])))
  print(paste("Number in:", nrow(hospitals[(hospitalid %in% hid)])))

  hospitals <- merge(hospitals, load_concepts("age", "eicu", verbose = FALSE), 
                     by = "patientunitstayid", all.x = T)


  res <- hospitals[hospitalid %in% hid & age >= 18L][["patientunitstayid"]]

  print(paste("Number of patients excluded on age:",
    nrow(hospitals[(hospitalid %in% hid) & (age < 18 | is.na(age))])))
  print(paste("Age exclusion from",
    length(unique(hospitals[(hospitalid %in% hid) & age < 18L][["hospitalid"]])), 
    "centres"))

  print(paste("Total Patients:", length(res)))

  res
}

cohort <- list(
  mimic = list(
    all = unique(id_col(mimic[age >= 18L])),
    bmi = unique(id_col(mimic[age >= 18L & !is.na(bmi)])),
    insulin = unique(id_col(mimic[age >= 18L & !is.na(bmi)]))
  ),
  eicu = list(
    all = pop_per_hosp(filter0),
    bmi = intersect(pop_per_hosp(filter0), 
                    unique(id_col(eicu[age >= 18L & !is.na(bmi)]))),
    insulin = intersect(pop_per_hosp(filter1), 
                        unique(id_col(eicu[age >= 18L & !is.na(bmi)])))
  ),
  hirid = list(
    all = unique(id_col(hirid[age >= 18L])),
    bmi = unique(id_col(hirid[age >= 18L & !is.na(bmi)])),
    insulin = unique(id_col(hirid[age >= 18L & !is.na(bmi)]))
  ),
  aumc = list(
    all = unique(id_col(aumc[age >= 18L])),
    bmi = unique(id_col(aumc[age >= 18L & !is.na(bmi)])),
    insulin = unique(id_col(aumc[age >= 18L & !is.na(bmi)]))
  )
)

config("cohort", cohort)

bmi_bins <- list(
  who = c(18.5, 25, 30, 35, 40),
  reg = seq(15, 40, 5)
)

config("bmi-bins", bmi_bins)

# generate the cohort for the Cox model
src <- c("mimic", "eicu", "hirid", "aumc")
cox_coh <- Map(function(x, y) setdiff(x[["insulin"]], remove_doi(y)), 
               cohort, src)
config("cox-cohort", cox_coh)









