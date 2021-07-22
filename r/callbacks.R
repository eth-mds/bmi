
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

DM_callback <- function(x, ...) {
  
  sub_var <- setdiff(names(x), meta_vars(x))
  if (sub_var == "icd9code") {
    
    x[, c(sub_var) := gsub(",.*", "", get(sub_var))]
    
  }
  
  intm <- data.frame(
    pid = id_col(x),
    icd9 = x[[setdiff(names(x), id_vars(x))]]
  )
  intm <- rowSums(comorbid_charlson(intm)[, c("DM", "DMcx")]) > 0
  
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
