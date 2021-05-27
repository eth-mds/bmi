library(ricu)
library(assertthat)
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
util_f <- file.path(root, "utils")
invisible(lapply(list.files(util_f), function(x) source(file.path(util_f, x))))

mmpp_mti <- function(source, concept, patient_ids = config("cohort")[[source]][["all"]]) {
  
  print("-------------------------")
  print(paste("source:", srcwrap(source), ", concept:", concept))
  
  xx <- load_concepts(concept, source, patient_ids = patient_ids, verbose = FALSE)
  xx <- merge(xx, stay_windows(source), all.x = T)
  xx <- xx[get(index_var(xx)) >= 0L & get(index_var(xx)) <= end]
  print(sprintf("Median number of measures per day of ICU stay %.2f", 
    median(xx[, .N / as.integer(end) * 24, by = eval(id_vars(xx))][["V1"]])))
  
  lact <- load_concepts("lact", source, patient_ids = patient_ids, verbose = FALSE)
  
  if(!is.null(patient_ids)) patient_ids <- intersect(
    patient_ids,
    unique(lact[[id_vars(lact)]])
  )
  
  x <- load_concepts(concept, source, patient_ids = patient_ids, verbose = FALSE)
  x <- merge(x, stay_windows(source), all.x = T)
  x <- x[get(index_var(x)) >= 0L & get(index_var(x)) <= end]
  
  mti <- function(t, e) {
    end <- unique(e)
    assertthat::assert_that(length(end) == 1L)
    t <- c(0, t, end)
    t_shift <- data.table::shift(t)
    #browser()
    median(t - t_shift, na.rm = T)
  }

  mti2 <- function(t, e) {
    end <- unique(e)
    assertthat::assert_that(length(end) == 1L)
    end / (length(t)+1)
  }
  
  intm <- x[, mti(get(index_var(x)), end), by = eval(id_vars(x))]
  print(length(unique(x[[id_vars(x)]])))
  plot(sort(x[, .N, by = eval(id_vars(x))][["N"]]), pch = 19, ylim = c(0, 20))
  print(paste("Median time interval between measurements:",
    median(x[, mti(get(index_var(x)), end), by = eval(id_vars(x))][["V1"]], na.rm = T))
  )
  
  print("-------------------------")

}

source <-  c("mimic", "eicu", "hirid")

sapply(source, mmpp_mti, "glu")
sapply(source, mmpp_mti, "lact")

hypo_summary <- function(source, patient_ids = config("cohort")[[source]][["all"]], ins_cf = 12L) {

  ins_hypo <- glycemia_treatment(source,
    vars = list(glu = list(time = 0L, imp_val = NA_real_), ins = list(time = ins_cf, imp_val = 0)),
    fill_na = T, patient_ids = patient_ids)

  ins_hypo <- ins_hypo[hypo == 1L]
  print("-------------------------")
  print(paste("Source:", srcwrap(source)))
  
  print(
    sprintf("%d episodes of hypoglycemia, %.2f%% of the cohort", 
      nrow(ins_hypo), 100 * nrow(ins_hypo) / length(patient_ids)
    )
  )

  print(paste("Median onset at", median(as.numeric(ins_hypo[[index_var(ins_hypo)]])), "hours"))
  print(paste("IQR of onset", quantile(as.numeric(ins_hypo[[index_var(ins_hypo)]]), c(0.25, 0.75)), "hours"))
  ppia <- 100* sum(ins_hypo[get(id_vars(ins_hypo)) %in% config("cohort")[[source]][["insulin"]]][["ins"]] > 0, na.rm = T) /
    nrow(ins_hypo[get(id_vars(ins_hypo)) %in% config("cohort")[[source]][["insulin"]]])
  print(paste("Proportion associated with insulin:", round(ppia, 2), "%"))
  print("-------------------------")
  return(T)
}

sapply(source, hypo_summary)

severe_hypo_summary <- function(source, patient_ids = config("cohort")[[source]][["all"]], ins_cf = 12L) {
  
  hypo_info <- glycemia_treatment(source,
    vars = list(glu = list(time = 0L, imp_val = NA_real_), ins = list(time = ins_cf, imp_val = 0)),
    fill_na = T, patient_ids = patient_ids, hypo.threshold = 2.2)
  
  ins_hypo <- hypo_info[hypo == 1L]
  
  assert_that(nrow(ins_hypo) == length(unique(id_col(ins_hypo))))
  
  print("-------------------------")
  
  print(paste("Source:", srcwrap(source)))
  
  print(
    sprintf("%d episodes of severe hypoglycemia, %.2f%%", 
      nrow(ins_hypo), 100 * nrow(ins_hypo) / length(patient_ids)
    )
  )
  
  print(paste("Median onset at", median(as.numeric(ins_hypo[[index_var(ins_hypo)]])), "hours"))
  
  print(paste("IQR of onset", quantile(as.numeric(ins_hypo[[index_var(ins_hypo)]]), c(0.25, 0.75)), "hours"))
  
  ppia <- 100* 
    sum(ins_hypo[get(id_vars(ins_hypo)) %in% config("cohort")[[source]][["insulin"]]][["ins"]] > 0, na.rm = T) / 
    nrow(ins_hypo[get(id_vars(ins_hypo)) %in% config("cohort")[[source]][["insulin"]]])
  
  print(sprintf("%.2f%% associated with insulin", ppia))
  
  hypo_info <- hypo_info[get(id_var(hypo_info)) %in% id_col(ins_hypo)]
  sbsh <- hypo_info[, min(c(head(glu, n = -1L), 108), na.rm = T), by = eval(id_var(hypo_info))][["V1"]]
  hpsh <- 100* sum(sbsh < 3.9 * 18.016) / nrow(ins_hypo)
  
  print(sprintf("In %.2f%% of cases severe hypo preceded by hypo", hpsh))
  
  print("-------------------------")
  return(T)
  
}

sapply(source, severe_hypo_summary)
