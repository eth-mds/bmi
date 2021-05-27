srcwrap <- function(src) {
  if (length(src) > 1) return(sapply(src, srcwrap))

  if(src == "mimic") {
    return("MIMIC-III")
  } else if (src == "eicu") {
    return("eICU")
  } else if (src == "hirid") {
    return("HiRID")
  } else if (src == "mimic_demo") {
    return("MIMIC-III (demo)")
  } else if (src == "eicu_demo") {
    return("eICU (demo)")
  } else if (src == "aumc") {
    return("AUMC")
  } else {
    return(src)
  }
}

undemo <- function(src) {

  if(src == "mimic_demo") return("mimic")
  if(src == "eicu_demo") return("eicu")

  return(src)

}


med_iqr <- function(x, patient_ids) {
  val_col <- setdiff(names(x), meta_vars(x))
  if(is_ts_tbl(x)) x <- x[get(index_var(x)) == 24L]
  quants <- quantile(x[[val_col]], probs = c(0.25, 0.5, 0.75), na.rm = T)
  res <- paste0(
    round(quants[2], 2), " (",
    round(quants[1], 2), "-",
    round(quants[3], 2), ")"
  )

  list(val_col, "Median (IQR)", res)
}

multi_med_iqr <- function(x, patient_ids) {

  val_cols <- setdiff(names(x), meta_vars(x))
  res <- lapply(
    val_cols, function(vcol) med_iqr(x[, c(meta_vars(x), vcol), with = FALSE], patient_ids)
  )

  lapply(1:3, function(i) {
    Reduce(c, lapply(res, `[[`, i))
  })

}

tab_design <- function(x, patient_ids) {

  val_col <- setdiff(names(x), meta_vars(x))
  res <- table(x[[val_col]])
  res <- round(100 * res / sum(res))

  if(val_col == "adm" & nrow(x) == 0L) {

    return(
      list(c("med", "surg", "other"), "%", rep(NA, 3))
    )

  }

  list(names(res), "%", as.integer(res))

}

percent_fun <- function(x, patient_ids) {

  val_col <- setdiff(names(x), meta_vars(x))

  if (val_col == "death") {

    return(list(val_col, "%", round(100 * sum(x[[val_col]]) / length(patient_ids))))

  }

  list(val_col, "%", round(100 * mean(x[[val_col]])))

}

concept_translator <- list(
  age = "Age (years)",
  med = "- Medical",
  other = "- Other",
  surg = "- Surgical",
  death = "Mortality",
  `Cohort size` = "Cohort size",
  los_icu = "ICU LOS",
  los_hosp = "Hospital LOS (days)",
  Male = "Gender (Male)",
  Female = "Gender (Female)",
  sofa = "- Total",
  sofa_resp_comp = "- Respiratory",
  sofa_coag_comp = "- Coagulation",
  sofa_cns_comp = "- CNS",
  sofa_liver_comp = "- Hepatic",
  sofa_cardio_comp = "- Cardiovascular",
  sofa_renal_comp = "- Renal"
)

charlson_callback <- function(x, ...) {
  
  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {
    
    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )

  intm <- rowSums(comorbid_charlson(intm))
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

elix_callback <- function(x, ...) {
  
  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {
    
    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  
  intm <- rowSums(comorbid_elix(intm))
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}
