
collect_hypo_cases <- function(tbl, max.hour = 240L) {
  # carry the hypo time backwards for 6 hours
  tbl <- slide(tbl, before = hours(0L), after = hours(5L), hypo_LA := max(hypo))

  # collect the relevant rows
  marks <- 6 * (1:(max.hour/6))
  rel_cols <- c(setdiff(names(tbl), meta_vars(tbl)))
  if (is.element("source", names(tbl))) rel_cols <- c(rel_cols, "source")
  collect <- tbl[get(index_var(tbl)) %in% marks, rel_cols, with = FALSE]
  collect[, hypo := NULL]
  collect <- data.table::setnames(collect, "hypo_LA", "hypo")

  collect

}

glycemia_treatment <- function(data_source,
  vars = list(glu = list(time = 24L, imp_val = NA_real_),
    lact = list(time = 24L, imp_val = 1), 
    ins_ifx = list(time = 12L, imp_val = 0),
    bmi = list(time = 0L, imp_val = NA), shock = list(time = 24L, imp_val = 0)), 
  fill_na = FALSE, patient_ids = NULL, hypo = TRUE, hypo.threshold = 3.9,
  sofa = FALSE, verbose = FALSE) {
  
  stat <- intersect(names(vars), c("bmi", "weight", "height", "age", "DM"))
  dyn <- setdiff(names(vars), stat)
  tbl <- load_concepts(dyn, data_source, patient_ids = patient_ids, 
                       verbose = verbose)
  static <- load_concepts(stat, data_source, patient_ids = patient_ids, 
                          verbose = verbose)

  # fill gaps
  tbl <- fill_gaps(tbl)
  tbl <- merge(tbl, static, all.x = TRUE)
  # reorder the columns appropriately
  tbl <- tbl[, c(meta_vars(tbl), names(vars)), with = FALSE]
  # carry-forward as specified
  tbl <- carry_values(tbl, vars)

  # fill NAs after the carry-forward if specified
  na_vals <- lapply(vars, function(x) x[["imp_val"]])
  if (fill_na) tbl <- tbl[, c(mget(meta_vars(tbl)),
                              Map(fill_missing, .SD, na_vals)),
                          .SDcols = names(vars)]
  
  if (sofa) {
    cmp <- "sofa_wo_cardio"
    sofa <- load_concepts(cmp, data_source, patient_ids = patient_ids, 
                          verbose = verbose)
    tbl <- merge(tbl, sofa, all.x = TRUE)
    if (fill_na) tbl <- replace_na(tbl, 0L, vars = cmp)
  }

  if (hypo) {
    # determine the first hypo time (glucose < 3.9 mmol/L)
    tbl[, hypo := as.integer(glu < hypo.threshold*18.016)]
    tbl[is.na(hypo), "hypo"] <- 0

    # delete everything after the first hypo time
    tbl <- tbl[, head(.SD, n = min(which(hypo == 1), length(hypo))),
               by = eval(id_vars(tbl))]
  }

  return(tbl)
}

tw_avg <- function(cnc, source, upto, hypo_censoring = TRUE,
                   patient_ids = config("cohort")[[source]][["bmi"]]) {

  x <- load_concepts(cnc, source, patient_ids = patient_ids, 
                     verbose = F)
  limits <- merge(x[, list(first_obs = min(get(index_var(x)))), 
                    by = c(id_vars(x))], 
                  stay_windows(source))
  limits[, start := pmin(start, first_obs)]
  if (hypo_censoring) {
    hg <- hypo(source, patient_ids, upto = hours(Inf))
    hg <- rename_cols(hg, "hypo_time", index_var(hg))
    hg <- as_id_tbl(hg)
    limits <- merge(limits, hg, all.x = TRUE)
    limits[is.na(hypo_time), hypo_time := hours(Inf)]
    limits[, end := pmin(end, hypo_time - hours(1L))]
  }
  
  x <- fill_gaps(x, limits = limits)
  x[, c(cnc) := data.table::nafill(get(cnc), "locf"), by = eval(id_vars(x))]
  
  x[get(index_var(x)) >= 0L & get(index_var(x)) <= upto]

  x[, list(target = mean(get(cnc), na.rm = T)), by = eval(id_vars(x))]

}

tw_avg_0imp <- function(cnc, source, upto, hypo_censoring = TRUE,
                       patient_ids = config("cohort")[[source]][["bmi"]]) {
  
  x <- load_concepts(cnc, source, patient_ids = patient_ids, 
                     verbose = F)
  limits <- merge(x[, list(first_obs = min(get(index_var(x)))), 
                    by = c(id_vars(x))], 
                  stay_windows(source))
  limits[, start := pmin(start, first_obs)]
  if (hypo_censoring) {
    hg <- hypo(source, patient_ids, upto = hours(Inf))
    hg <- rename_cols(hg, "hypo_time", index_var(hg))
    hg <- as_id_tbl(hg)
    limits <- merge(limits, hg, all.x = TRUE)
    limits[is.na(hypo_time), hypo_time := hours(Inf)]
    limits[, end := pmin(end, hypo_time - hours(1L))]
  }
  
  x <- fill_gaps(x, limits = limits)
  x[is.na(get(cnc)), c(cnc) := 0]
  
  x[get(index_var(x)) >= 0L & get(index_var(x)) <= upto]
  
  targ <- x[, list(target = mean(get(cnc), na.rm = T)), by = eval(id_vars(x))]
  
  ful <- id_tbl(patient_ids)
  ful <- rename_cols(ful, id_var(targ), id_var(ful))
  ful <- merge(ful, targ, all.x = TRUE)
  ful[is.na(target), target := 0]
  
  ful
  
}

med_dur <- function(cnc, source, upto, hypo_censoring = TRUE,
                   patient_ids = config("cohort")[[source]][["bmi"]]) {
  
  x <- load_concepts(cnc, source, patient_ids = patient_ids, 
                     verbose = F)
  limits <- merge(x[, list(first_obs = min(get(index_var(x)))), 
                    by = c(id_vars(x))], 
                  stay_windows(source))
  limits[, start := pmin(start, first_obs)]
  if (hypo_censoring) {
    hg <- hypo(source, patient_ids, upto = hours(Inf))
    hg <- rename_cols(hg, "hypo_time", index_var(hg))
    hg <- as_id_tbl(hg)
    limits <- merge(limits, hg, all.x = TRUE)
    limits[is.na(hypo_time), hypo_time := hours(Inf)]
    limits[, end := pmin(end, hypo_time - hours(1L))]
  }
  
  if (is_win_tbl(x)) x <- expand(x, aggregate = TRUE)
  if (is.logical(x[[cnc]])) x[, c(cnc) := as.double(get(cnc))]
  x <- fill_gaps(x, limits = limits)
  x[is.na(get(cnc)), c(cnc) := 0]
  
  x[get(index_var(x)) >= 0L & get(index_var(x)) <= upto]
  
  targ <- x[, list(target = mean(get(cnc), na.rm = T)), by = eval(id_vars(x))]
  
  ful <- id_tbl(patient_ids)
  ful <- rename_cols(ful, id_var(targ), id_var(ful))
  ful <- merge(ful, targ, all.x = TRUE)
  ful[is.na(target), target := 0]
  
  ful
  
}

glu_freq <- function(source, upto,
                     patient_ids = config("cohort")[[source]][["bmi"]]) {

  gluc <- load_concepts("glu", source, patient_ids = patient_ids)
  gluc <- gluc[index_col(gluc) >= 0L & index_col(gluc) <= upto]

  gluc <- gluc[, list(n_measures = .N), by = eval(id_var(gluc))]

  gluc <- merge(gluc, load_concepts("los_icu", source, 
                                    patient_ids = patient_ids))

  gluc[, los_icu := min(los_icu, upto / 24), by = eval(id_var(gluc))]
  gluc[, target := los_icu * 24 / n_measures]

  gluc[, c(id_var(gluc), "target"), with = F]

}

get_target <- function(source, target, upto, patient_ids) {

  if (target == "tw_avg_glucose") {

    res <- tw_avg("glu", source, upto, patient_ids = patient_ids)

  } else if (target == "dur_TPN") {
    
    res <- med_dur("TPN", source, upto, patient_ids = patient_ids)
    
  } else if (target == "dur_enteral") {
    
    res <- med_dur("enteral", source, upto, patient_ids = patient_ids)
    
  } else if (target == "dur_cortico") {
    
    res <- med_dur("cortico", source, upto, patient_ids = patient_ids)
    
  } else if (target == "tw_avg_dextrose") {
    
    res <- tw_avg_0imp("dex_amount", source, upto, patient_ids = patient_ids)
    
  } else if (target == "min_pafi") {
    
    res <- w_value(source, "pafi", dir = "decreasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "max_insulin") {

    res <- w_value(source, "ins_ifx", dir = "increasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")

  } else if (target == "max_crea") {
    
    res <- w_value(source, "crea", dir = "increasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "max_lact") {
    
    res <- w_value(source, "lact", dir = "increasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "max_ast") {
    
    res <- w_value(source, "ast", dir = "increasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "min_plt") {
    
    res <- w_value(source, "plt", dir = "decreasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "min_pafi") {
    
    res <- w_value(source, "pafi", dir = "decreasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "min_map200") {
    
    res <- w_value(source, "map_beta200", dir = "decreasing", upto = upto,
                   patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    
  } else if (target == "max_insulin_wnorm") {

    res <- w_value(source, "ins_ifx", dir = "increasing", upto = upto,
                   patient_ids = patient_ids)
    res <- merge(res, load_concepts("weight", source), patient_ids = patient_ids)
    res <- rename_cols(res, "target", "w_val")
    res[, target := target/weight]

  } else if (target == "glu_freq") {

    res <- glu_freq(source, upto, patient_ids)

  } else if (target == "fhm") {

    res <- hypo(source, patient_ids, upto = upto, value = T)
    res <- rename_cols(res, "target", "glu")

  } else if (target == "hypo") {
    
    res <- hypo(source, patient_ids, upto = upto)
    res <- rename_cols(res, "hypo", "hg")
    
  } else res <- load_concepts(target, source, patient_ids = patient_ids)

  if (is.element("target", names(res))) res <- rename_cols(res, target, "target")

  res

}

high_freq <- function(src, patient_ids, upto) {
  
  tbl <- get_target(src, "glu_freq", upto = hours(Inf), 
                    patient_ids = patient_ids)
  fthres <- quantile(tbl[["glu_freq"]], 0.25)
  id_col(tbl[glu_freq < fthres])
  
}

hypo_only <- function(src, patient_ids, upto) {
  
  id_col(hypo(src, patient_ids, upto = upto))
  
}

assoc_hmap <- function(res, gm = "gv_cv", bins = c(10, 20, 30),
                       bin_lvls = c("<10%", "10-20%", "20-30%", ">30%"),
                       y_label = "Glucose coefficient of variation",
                       title = "Glucose variability and mortality") {
  
  bern_ci <- function(p, n) {
    if (length(p) > 1) return(unlist(Map(bern_ci, p, n)))
    
    if (n > 100) {
      sgm <- sqrt(p * (1-p) / n)
      if (p - 1.96*sgm >= 0) {
        return(paste0("(", spec_dec(p - 1.96*sgm, 2), "-", 
                      spec_dec(p + 1.96*sgm, 2), ")"))
      }
    }
    pqs <- quantile(rbinom(100, n, p)/n, c(0.025, 0.975))
    return(paste0("(", spec_dec(pqs[1], 2), "-", spec_dec(pqs[2], 2), ")"))
  }
  
  res <- res[!is.na(get(gm))]
  res[, 
      c(gm) := factor(.bincode(get(gm), c(-Inf, bins, Inf)), 
                      labels = bin_lvls)
  ]
  res[, DM := factor(DM, labels = c("No DM", "DM"))]
  
  dat <- res[, list(Mortality = mean(death), gsize = .N), 
             by = c(gm, "bmi_bins", "DM")]
  dat[, 
      txt_lab := paste0(spec_dec(Mortality, 2), "\n", bern_ci(Mortality, gsize))
     ]
  
  ggplot(dat, aes_string(x = "DM", y = gm, fill = "Mortality")) +
    geom_bin_2d(aes_string(y = gm)) +
    facet_grid(cols = vars(bmi_bins)) + 
    scale_fill_viridis_c(limits = c(0.00, 0.26), breaks = c(0.05, 0.15, 0.25)) +
    geom_text(aes(label = txt_lab), color = "red", size = 3) +
    xlab("BMI bins") + ylab(y_label) +
    theme_bw() + ggtitle(title) + xlab(NULL) +
    theme(legend.position = "bottom")
  
}
