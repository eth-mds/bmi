
liver_damage <- function(data_source = "mimic", id_type = "icustay",
                         patient_ids = NULL, verbose = FALSE){

  liv_dam <- load_concepts(
    c("bili", "alt", "ast"), src = data_source,
    id_type = id_type, patient_ids = patient_ids, verbose = verbose
  )
  liv_dam[, liver_damage := as.integer(bili > 2 | alt > 45 | ast > 45)]
  liv_dam <- liv_dam[, c(meta_vars(liv_dam), "liver_damage"), with = FALSE]

  liv_dam
}

shock <- function(data_source = "mimic", id_type = "icustay",
                  patient_ids = NULL, verbose = FALSE) {

  x <- load_concepts(c("dopa_rate", "norepi_rate", "dobu_rate", "epi_rate", "map"),
    data_source, id_type = id_type, patient_ids = patient_ids, verbose = verbose)

  if (data_source == "hirid") x[, dopa := NA]

  x <- x[, c(meta_vars(x), "dopa_rate", "norepi_rate", "dobu_rate", "epi_rate",
             "map"), with = FALSE]
  x[is.na(map), "map"] <- 100

  shock <- x[, as.integer(any(!is.na(dopa_rate), !is.na(dobu_rate),
                              !is.na(epi_rate), !is.na(norepi_rate), map < 60)),
             by = eval(meta_vars(x))]
  shock <- data.table::setnames(shock, "V1", "shock")
  shock <- shock[shock == 1]

  return(shock)
}

mech_vent <- function(data_source = "mimic", id_type = "icustay",
                      patient_ids = NULL, verbose = FALSE) {

  mechv <- ricu:::sofa_vent(
    load_concepts("vent_start", data_source, id_type = id_type,
                  patient_ids = patient_ids, verbose = verbose),
    load_concepts("vent_end", data_source, id_type = id_type,
                  patient_ids = patient_ids, verbose = verbose),
    hours(6L), hours(2L), hours(1L)
  )

  mechv[, mech_vent := as.integer(vent)]
  mechv <- mechv[, c(meta_vars(mechv), "mech_vent"), with = FALSE]

  return(mechv)
}

load_ood <- function(data_source, concepts, id_type = "icustay", patient_ids = NULL) {

  res <- list()

  for(i in 1:length(concepts)) {

    f <- eval(parse(text = (concepts[i])))
    res[[i]] <- f(data_source, id_type = id_type, patient_ids = patient_ids)

  }

  if (length(res) == 1) return(res[[1]])

  res <- Reduce(function(x, y) merge(x, y, all = TRUE), res)

  return(res)

}

collect_hypo_cases <- function(tbl, max.hour = 240L) {
  # carry the hypo time backwards for 6 hours
  tbl <- slide(tbl, before = hours(0L), after = hours(5L), hypo_LA := max(hypo))

  # collect the relevant rows
  marks <- 6 * (1:(max.hour/6))
  rel_cols <- c(setdiff(names(tbl), meta_vars(tbl)))
  collect <- NULL
  for (mark in marks) {
    tmp <- tbl[get(index_var(tbl)) == mark]
    collect <- rbind(
      collect,
      as.matrix(tmp[, rel_cols, with = FALSE])
    )
  }

  collect <- data.table::data.table(collect)
  collect[, hypo := NULL]
  collect <- data.table::setnames(collect, "hypo_LA", "hypo")

  collect

}

glycemia_treatment <- function(data_source,
  vars = list(glu = list(time = 24L, imp_val = NA_real_),
    lact = list(time = 24L, imp_val = 1), ins = list(time = 12L, imp_val = 0),
    shock = list(time = 24L, imp_val = 0),
    liver_damage = list(time = 48L, imp_val = 0)), fill_na = FALSE,
    id_type = "icustay", patient_ids = NULL, hypo = TRUE, hypo.threshold = 3.9,
    sofa = TRUE, verbose = FALSE) {

  dict <- ricu::get_config("concept-dict", ricu:::default_config_path())

  in_dict <- intersect(names(vars), names(dict))
  out_dict <- setdiff(names(vars), in_dict)
  tbl1 <- load_concepts(in_dict, data_source, id_type = id_type,
                        patient_ids = patient_ids, verbose = verbose)

  if (length(out_dict) > 0L) {
    tbl2 <- load_ood(data_source, out_dict, id_type = id_type,
                     patient_ids = patient_ids)
    tbl <- merge(tbl1, tbl2, all = TRUE)
  } else {
    tbl <- tbl1
  }

  if (data_source == "mimic" & is.element("ins", names(tbl))) {
    tbl[ins == 0, "ins"] <- 2 # MIMIC carevue imputation
  } 
  # reorder the columns appropriately
  tbl <- tbl[, c(meta_vars(tbl), names(vars)), with = FALSE]

  # fill gaps
  tbl <- fill_gaps(tbl)
  tbl <- carry_values(tbl, vars)

  # fill NAs after the carry-forward if specified
  na_vals <- lapply(vars, function(x) x[["imp_val"]])
  if (fill_na) tbl <- tbl[, c(mget(meta_vars(tbl)),
                              Map(fill_missing, .SD, na_vals)),
                          .SDcols = names(vars)]
  
  if (sofa) {
    sofa <- load_concepts("sofa", data_source, keep_components = TRUE,
                          verbose = verbose)
    sofa[, c("sofa_liver_comp", "sofa_cardio_comp") := NULL]
    sofa <- replace_na(sofa, 0L)
    
    cmp <- c("sofa_coag_comp", "sofa_renal_comp", "sofa_cns_comp", 
             "sofa_resp_comp")
    # cmp <- "sofa_mlc"
    # sofa[, sofa_mlc := sofa_coag_comp + sofa_renal_comp + sofa_cns_comp +
    #                    sofa_resp_comp]
    sofa <- sofa[, c(meta_vars(sofa), cmp), with = FALSE]
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

tw_avg_glucose <- function(source, upto,
                           patient_ids = config("cohort")[[source]][["bmi"]]) {

  x <- fill_gaps(load_concepts("glu", source, patient_ids = patient_ids, 
                               verbose = F))
  x[, glu := data.table::nafill(glu, "locf"), by = eval(id_vars(x))]

  wins <- stay_windows(source)
  x <- merge(x, wins)
  x[get(index_var(x)) >= 0L & get(index_var(x)) <= upto]

  x[, list(target = mean(glu, na.rm = T)), by = eval(id_vars(x))]

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

    res <- tw_avg_glucose(source, upto, patient_ids)

  } else if (target == "max_insulin") {

    res <- w_value(source, "ins", dir = "increasing", upto = upto,
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

    res <- w_value(source, "ins", dir = "increasing", upto = upto,
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

cox_treatment <- function(data_source,
                          vars = list(
                            glu = list(type = "locf", imp_val = NA,
                                       before_val = 108),
                            lact = list(type = "locf", imp_val = NA,
                                        before_val = 1),
                            ins = list(type = "const", imp_val = 0,
                                       before_val = 0),
                            shock = list(type = "locf", imp_val = NA,
                                         before_val = 0),
                            liver_damage = list(type = "locf", imp_val = NA,
                                                before_val = 0L),
                            death = list(type = "const", imp_val = FALSE,
                                         before_val = 0L)),
                          id_type = "icustay", patient_ids = NULL, fill_gaps = T,
                          hypo.threshold = 3.9, verbose = FALSE) {

  dict <- ricu::get_config("concept-dict", ricu:::default_config_path())

  in_dict <- intersect(names(vars), names(dict))
  out_dict <- setdiff(names(vars), in_dict)
  tbl1 <- load_concepts(in_dict, data_source, id_type = id_type,
                        patient_ids = patient_ids, verbose = verbose)

  if (length(out_dict) > 0L) {
    tbl2 <- load_ood(data_source, out_dict, id_type = id_type,
                     patient_ids = patient_ids)
    tbl <- merge(tbl1, tbl2, all = TRUE)
  } else {
    tbl <- tbl1
  }

  if (data_source == "mimic" & is.element("ins", names(tbl))) tbl[ins == 0,
                                                                  "ins"] <- 2
  # reorder the columns appropriately
  tbl <- tbl[, c(meta_vars(tbl), names(vars)), with = FALSE]

  # fill gaps
  if (fill_gaps) tbl <- fill_gaps(tbl)
  dat <- replace_na(tbl, sapply(vars, `[[`, "imp_val"),
                    type = sapply(vars, `[[`, "type"),
                    by_ref = TRUE,
                    vars = names(vars),
                    by = id_vars(tbl))


  dat <- replace_na(dat, sapply(vars, `[[`, "before_val"),
                    type = rep("const", length(vars)),
                    by_ref = TRUE,
                    vars = names(vars),
                    by = id_vars(dat))

  dat <- merge(dat, load_concepts("bmi", data_source, verbose = verbose),
               all.x = T)

  dat <- merge(dat, stay_windows(data_source), all.x = T)
  dat <- dat[get(index_var(dat)) >= start &
             get(index_var(dat)) <= end]
  dat <- dat[, c("start", "end") := NULL]

  ns_coh <- unique(id_col(dat[death == T]))
  dat <- dat[,
             end := shift(get(index_var(dat)), n = -1L, fill = NA,
                            type = "lag"), by = c(id_vars(dat))]

  dat <- dat[get(id_vars(dat)) %in% ns_coh,
             death := shift(death, n = -1L, fill = NA,
                            type = "lag"), by = c(id_vars(dat))]


  # fill NAs after locf

  dat <- dat[, head(.SD, n = match(TRUE, death, .N)), by = c(id_vars(dat))]

  # get hypo
  dat[, hypo := as.integer(glu <= hypo.threshold * 18.016)]
  dat[, hypo := cummax(hypo), by = c(id_vars(dat))]

  # scale glucose for better interpretability
  dat <- dat[, glu := glu / 18.016]

  return(dat)
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
