
ts_to_win_2hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(120L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

ts_to_win_6hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

mimic_presc_cort <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(index_var(x)) := get(index_var(x)) + hours(9L)]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

hirid_pharma_win6 <- function(x, dur_var, group_var, ...) {
  
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  
  x[, c(dur_var) := max(get(index_var(x))) - min(get(index_var(x))) + hours(6L), 
    by = c(group_var)]
  x <- x[, head(.SD, n = 1L), by = c(group_var)]
  x[, c(dur_var) := `units<-`(get(dur_var), "mins")]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

hirid_pharma_win2 <- function(x, dur_var, group_var, ...) {
  
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  
  x[, c(dur_var) := max(get(index_var(x))) - min(get(index_var(x))) + hours(2L), 
    by = c(group_var)]
  x <- x[, head(.SD, n = 1L), by = c(group_var)]
  x[, c(dur_var) := `units<-`(get(dur_var), "mins")]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

aumc_cortico <- function(x, dur_var, ...) {
  
  x[, c(dur_var) := get(dur_var) + mins(360L)]
  
}

DM910_callback <- function(x, val_var, ...) {
  
  if (val_var == "icd9code") {
    
    x[, c(val_var) := gsub(",.*", "", get(val_var))]
    
  }
  
  DM_map <- list(
    DM = c(icd9_map_charlson$DM, icd10_map_charlson$DM),
    DMcx = c(icd9_map_charlson$DMcx, icd10_map_charlson$DMcx)
  )
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(comorbid(intm, map = DM_map)[, c("DM", "DMcx")]) > 0
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

SMK_callback <- function(x, val_var, ...) {
  
  if (val_var == "icd9code") {
    
    x[, c(val_var) := gsub(",.*", "", get(val_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(
    icd9_comorbid(intm, map = list(Smoking = c("3051", "305.1")))[, c("Smoking"), 
                                                      drop = FALSE]) > 0
  
  res <- id_tbl(
    id = as.integer(names(intm)),
    val = intm, id_vars = "id"
  )
  
  names(res) <- names(x)
  
  res
}

map_beta_200 <- function (..., match_win = hours(2L), beta = 200, 
                          interval = NULL) {
  
  cnc <- c("map", "norepi_equiv")
  res <- ricu:::collect_dots(cnc, interval, ...)
  
  assert_that(ricu:::is_interval(match_win), 
              match_win > ricu:::check_interval(res))
  
  on12 <- paste(meta_vars(res[[1L]]), "==", meta_vars(res[[2L]]))
  on21 <- paste(meta_vars(res[[2L]]), "==", meta_vars(res[[1L]]))
  res <- rbind(res[[1L]][res[[2L]], on = on12, roll = match_win],
               res[[2L]][res[[1L]], on = on21, roll = match_win])
  res <- unique(res)
  res[is.na(get(cnc[2L])), cnc[2L]] <- 0 # impute 0 for vasos where needed
  
  res <- res[!is.na(get(cnc[1L])) & !is.na(get(cnc[2L])), ]
  res <- res[, `:=`(c(paste0("map_beta", beta)), get(cnc[1L]) - 
                      beta*get(cnc[2L]))]
  res <- rm_cols(res, cnc)
  res
}

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

dex_amount_callback <- function(...) {
  x <- list(...)[["dex"]]
  ivl <- list(...)[["interval"]]
  c_fct <- as.double(ivl, units = units(x$dur_var))
  # make into amount
  x[, dex := dex * as.double(dur_var) / c_fct]
  x[, distr := as.double(ricu:::re_time(dur_var, ivl)) + 1]
  x[, dex := dex / distr]
  x[, distr := NULL]
  x <- rename_cols(x, "dex_amount", "dex")
  
  expand(x, aggregate = "sum")
}

liver_damage_callback <- function(..., interval) {
  
  liv_dam <- Reduce(function(x, y) merge(x, y, all = TRUE), list(...))
  liv_dam[, liver_damage := as.integer(bili > 2 | alt > 45 | ast > 45)]
  liv_dam[is.na(liver_damage), liver_damage := 0]
  liv_dam <- liv_dam[, c(meta_vars(liv_dam), "liver_damage"), with = FALSE]
  
  liv_dam
  
}

shock_callback <- function(..., interval) {
  
  shock <- Reduce(function(x, y) merge(x, y, all = TRUE), list(...))
  shock[, shock := as.integer(any(!is.na(dopa_rate), !is.na(dobu_rate),
                                  !is.na(epi_rate), !is.na(norepi_rate), 
                                  map < 60)),
        by = eval(meta_vars(shock))]
  shock[, c(meta_vars(shock), "shock"), with = FALSE]
  
}

mimic_version <- function(x, val_var, ...) {
  
  x[, c(val_var) := ifelse(get(val_var) == "metavision", "new", "old")]
  
}

aumc_version <- function(x, val_var, ...) {
  
  x[, c(val_var) := ifelse(get(val_var) == "2010-2016", "new", "old")]
  
}

gen_version <- function(x, val_var, ...) {
  
  x[, c(val_var) := "new"]
  
}

ins_ifx_cb <- function(ins, interval) {
  
  if (length(id_vars(ins)) == 1L & id_vars(ins)[1L] == "icustay_id") {
    ins[ins == 0, ins := 2]
  }
  
  if (length(id_vars(ins)) > 1L) {
    ins[ins == 0 & grepl("mimic", source), ins := 2]
  } 
  
  rename_cols(ins, "ins_ifx", "ins")
  
}

sofa_woc <- function (..., worst_val_fun = max_or_na, explicit_wins = FALSE, 
                      win_length = hours(24L), keep_components = FALSE, 
                      interval = NULL) 
{
  cnc <- c("sofa_resp", "sofa_coag", "sofa_liver", "sofa_cns", "sofa_renal")
  dat <- ricu:::collect_dots(cnc, interval, ..., merge_dat = TRUE)
  expr <- substitute(lapply(.SD, fun), list(fun = worst_val_fun))
  if (isFALSE(explicit_wins)) {
    res <- fill_gaps(dat)
    res <- slide(res, !!expr, before = win_length, full_window = FALSE, 
                 .SDcols = cnc)
  }
  else {
    if (isTRUE(explicit_wins)) {
      assert_that(ricu:::is_scalar(win_length), ricu:::is_interval(win_length))
      ind <- index_var(dat)
      win <- dat[, list(max_time = max(get(ind))), by = c(id_vars(dat))]
      win <- win[, `:=`(c("min_time"), get("max_time") - 
                          win_length)]
      res <- hop(dat, !!expr, win, .SDcols = cnc)
    }
    else {
      res <- slide_index(dat, !!expr, explicit_wins, before = win_length, 
                         full_window = FALSE, .SDcols = cnc)
    }
  }
  res <- res[, `:=`(c("sofa_wo_cardio"), rowSums(.SD, na.rm = TRUE)), 
             .SDcols = cnc]
  if (isTRUE(keep_components)) {
    res <- rename_cols(res, paste0(cnc, "_comp"), cnc, by_ref = TRUE)
  }
  else {
    res <- rm_cols(res, cnc, by_ref = TRUE)
  }
  res
}

is_hypo_cb <- function(glu, interval, ...) {
  
  glu[, list(is_hypo = any(glu <= 70)), by = c(id_vars(glu))]
  
}

bin_bmi <- function(bmi, ...) {
  
  breaks <- c(-Inf, config("bmi-bins")[["who"]], Inf)
  
  bmi[, bmi_bins := factor(.bincode(bmi, breaks))]
  levels(bmi[["bmi_bins"]]) <- bin_labels(config("bmi-bins")[["who"]], "kg/m2")
  
  id_var <- id_vars(bmi)
  bmi[, c(id_var, "bmi_bins"), with = F]
  
}

gv_cv <- function(glu, ...) {
  ind <- index_var(glu)
  glu <- glu[get(ind) >= hours(0L)]
  
  glu[, list(gv_cv = 100 * sd(glu) / mean(glu)), by = c(id_vars(glu))]
}

gv_sd <- function(glu, ...) {
  ind <- index_var(glu)
  glu <- glu[get(ind) >= hours(0L)]
  
  glu[, list(gv_sd = sd(glu)), by = c(id_vars(glu))]
}

hypo_term <- function(x, min_dur, max_dur) {
  
  mrg <- function(x) cumsum(c(TRUE, x[-length(x)]))
  
  trm <- function(x, mx) {
    ind <- max(which(x > 0L))
    sft <- ind + mx
    if (sft <= length(x)) replace(x, sft, -x[ind]) else x
  }
  
  idv <- id_vars(x)
  idx <- index_var(x)
  
  assert_that(!any(c("diff", "hdif", "nhyp", "impu") %in% colnames(x)))
  
  x <- x[!is.na(hypo), diff := c(diff(get(idx)), Inf), by = c(idv)]
  x <- x[is_true(hypo > 0L), hdif := c(diff(get(idx)), Inf), by = c(idv)]
  x <- x[!is.na(hdif), nhyp := .N, by = c(idv)]
  
  x <- x[is_true(nhyp > 1L), hypo := mrg(hdif > min_dur | diff < hdif),
         by = c(idv)]
  x <- x[, c("diff", "hdif", "nhyp") := NULL]
  
  x <- x[, impu := data.table::nafill(hypo, "locf"), by = c(idv)]
  x <- x[, impu := data.table::nafill(impu, fill = 0)]
  
  mxd <- as.integer(
    ceiling(max_dur / as.double(interval(x), units = units(max_dur)))
  )
  
  x <- x[impu > 0L, hypo := trm(hypo, mxd), by = c(idv, "impu")]
  x <- x[, impu := NULL]
  
  x
}

#' @return Constructed from a `hypo` column a returned by `hypo_term()`,
#' several columns are added to the passed `ts_tbl` and returned as such:
#' * `hypo_imp`: using an locf imputation scheme, hypo periods are marked by
#' ascending even numbers and the preceding non-hypo periods by odd integers.
#' * `hypo_epi`: hypo episodes are constructed from the `hypo_imp` column such
#' that a given hypo periods and the stretch leading up to it have assigned the
#' same integer.
#' * `start_time`: Time relative to the start of a given hypo episode.
#' * `onset_time`: Either `NA` if a given hypo episode does not contain a hypo
#' onset or the time relative to hypo onset.
hypo_augm <- function(x) {
  
  shift <- function(x) data.table::fifelse(x > 0L, x * 2L, x * -2L + 1L, 1L)
  
  timed <- function(tim, i, j) list(tim - tim[i], tim - tim[j])
  
  idv <- id_vars(x)
  
  assert_that(!"temp" %in% colnames(x))
  
  x <- x[, c("temp", "hypo") := list(
    data.table::fifelse(hypo < 0L, 0L, hypo), NULL
  )]
  
  x <- x[, temp := data.table::nafill(temp, "locf"), by = c(idv)]
  x <- x[, temp := data.table::nafill(temp, fill = 0L)]
  
  x <- x[, temp := c(0L, diff(temp)), by = c(idv)]
  x <- x[, temp := data.table::fifelse(
    temp == 0L, NA_integer_, temp
  )]
  x <- x[, temp := data.table::nafill(temp, "locf"), by = c(idv)]
  x <- x[, c("hypo_imp", "temp") := list(shift(temp), NULL)]
  x <- x[, hypo_epi := (hypo_imp + 1L) %/% 2L]
  
  x <- x[, c("start_time", "onset_time") := timed(
    get(index_var(x)), 1L, c(FALSE, diff(hypo_imp) > 0L)),
    by = c(idv, "hypo_epi")
  ]
  
  x
}

hypo_cb <- function(glu, ...) {
  
  onset_id <- function(x) replace(x, x, seq_len(sum(x)))
  
  glu_var <- data_var(glu)
  id_vars <- id_vars(glu)
  
  glu <- glu[, c("tmp") := is_true(get(glu_var) <= 3.9 * 18.016)]
  glu <- glu[, c("hypo") := onset_id(get("tmp")), by = c(id_vars)]
  glu <- glu[, c(glu_var, "tmp") := NULL]
  glu <- fill_gaps(glu)
  
  glu
}

hypo_episode <- function(hypo, min_dur = hours(6L), max_dur = hours(4L), ...) {
  
  hypo <- hypo_term(hypo, min_dur, max_dur)
  hypo <- hypo_augm(hypo)
  
  hypo[, c("hypo_epi", "start_time", "onset_time") := NULL]
  
  rename_cols(hypo, "hypo_epi", "hypo_imp")
}

hypo_cnt <- function(hypo_epi, ...) {
  hypo_epi[, list(hypo_cnt = floor(max(hypo_epi / 2))), by = c(id_vars(hypo_epi))]
}

hypo_dur <- function(glu, ... ) {
  
  ind <- index_var(glu)
  glu <- glu[get(ind) >= hours(0L)]
  glu[, list(hypo_dur = mean(glu <= 70)), by = c(id_var(glu))]
  
}

hypo_sev <- function(hypo_epi, glu, ...) {
  
  hypo_epi <- merge(hypo_epi, glu, all.x = TRUE)
  
  hypo_epi[hypo_epi %% 2 == 0, 
           list(hypo_sev = min(glu, na.rm = TRUE)), 
           by = c(id_var(hypo_epi), "hypo_epi")]
  
}

tw_avg_gluc <- function(glu, icu_end, upto = hours(Inf), hypo_censoring = TRUE, ...) {

  ind <- index_var(glu)
  limits <- merge(glu[, list(first_obs = min(get(ind))), 
                    by = c(id_vars(glu))], 
                  icu_end)
  limits[, start := pmin(hours(0L), first_obs)]
  if (hypo_censoring) {
    hg <- glu[, list(hypo_time = head(get(ind)[glu <= 70], 1L)), 
              by = c(id_vars(glu))]
    limits <- merge(limits, hg, all.x = TRUE)
    limits[is.na(hypo_time), hypo_time := hours(Inf)]
    limits[, end := pmin(icu_end, hypo_time - hours(1L))]
  }
  
  glu <- fill_gaps(glu, limits = limits)
  glu[, glu := data.table::nafill(glu, "locf"), by = eval(id_vars(glu))]
  
  glu[get(ind) >= 0L & get(ind) <= upto]
  
  glu[, list(tw_avg_glu = mean(glu, na.rm = T)), by = eval(id_vars(glu))]
  
}

stay_win_cb <-function (x, id_type, interval)
{
  assert_that(id_type == "icustay") # assumes only icustay
  cfg <- as_id_cfg(x)
  res <- id_map(x, id_vars(cfg[id_type]), id_vars(cfg[id_type]), 
                NULL, "end")
  res <- res[, `:=`(c("val_var", "end"), list(get("end"), NULL))]
  res <- change_interval(res, interval)
}
