
bin_labels <- function(breaks, unit, lower0 = TRUE) {

  x_labels <- sapply(1:(length(breaks)-1),
    function(x) paste0("[", breaks[x], "-", breaks[x+1], "]")
  )
  first_label <- paste0("< ", breaks[1])
  if (lower0) first_label <- paste0("[0-", breaks[1], "]")
  x_labels <- c(first_label, x_labels, paste0("> ", breaks[length(breaks)]))
  x_labels <- paste(x_labels, unit)

  return(x_labels)
}

theme_fp <- function(...) {
  theme_bw(...) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

reshape_frame <- function(frame, break_list) {

  for(col in names(break_list)) {

    breaks <- break_list[[col]]

    x <- frame[, col]
    x <- as.factor(.bincode(x, breaks = c(-Inf, sort(breaks), Inf)))
    if(length(breaks) == 1 & breaks[1] == 0.5) {
      levels(x) <- c("_no", "_yes")
    } else {
      levels(x) <- c(
        paste0(" [", c(0, breaks[-length(breaks)]), "-", breaks, "]"),
        paste(" >", breaks[length(breaks)])
      )
    }
    frame[, col] <- x

  }

  return(frame)

}

extract_score <- function(model, p.threshold = 0.05) {

  score <- model$coefficients * as.integer(summary(logit)$coefficients[, 4] < p.threshold)
  score <- score[-1]
  score <- score / min(abs(score[score != 0]))
  score <- round(score)

  score

}

discretize <- function(x, breaks) {

  if (is.null(breaks)) return(x)
  return(.bincode(x, c(-Inf, breaks, Inf)))

}

fill_missing <- function(x, val) {
  x[is.na(x)] <- val
  return(x)
}

carry_fwd <- function(x, time = 1L) {

  res <- x[length(x)]

  if(is.na(res)) {

    cand <- which(!is.na(x))
    if(length(cand) == 0) return(res)
    cand <- max(cand)
    if ((length(x) - cand) <= time) res <- x[cand]

  }

  res

}

carry_bwd <- function(x, time = 1L) {

  res <- x[1]


  if (is.na(res)) {
    cand <- which(!is.na(x))
    if(length(cand) == 0) return(res)
    cand <- min(cand)
    if (cand <= time) res <- x[cand]
  }

  res

}

median_na <- function(x) {
  if(all(is.na(x))) return(x[1])
  return(median(x, na.rm = T))
}

make_plots <- function(res, save_plots, folder, plot_name, bottom = NULL,
                       width = 12, height = 5, dpi.res = 300, n.cols = length(res)) {

  wd <- getwd()

  ranges <- lapply(res, function(x) ggplot_build(x)$layout$panel_scales_y[[1]]$range$range)
  y_limits <- Reduce(function(x, y) c(min(x, y), max(x, y)), ranges)
  res <- lapply(res, function(x) x + ylim(y_limits))
  plot <- plot_grid(plotlist = res,
    ncol = n.cols, labels = c("A", "B", "C", "D")[1:length(res)])
  grid.arrange(arrangeGrob(plot, bottom = bottom))

  if(save_plots) {

    ggsave(file.path(wd, "paper", "figures", folder, paste0(plot_name, ".tif")),
      device = "tiff", width = width, height = height, dpi = dpi.res)

  }

}

carry_values_single <- function(tbl, col, time, dir = "forward") {
  assert_that(has_no_gaps(tbl))

  if (time == 0L) return(tbl)

  tbl <- data.table::copy(tbl)
  on.exit(tbl[, lag := NULL])

  if(is.null(dir)) dir <- "forward"
  dir <- as.integer((-1)^(dir == "backward"))
  for(i in seq.int(1, time)) {
    tbl[, lag := data.table::shift(.SD, n = dir), by = eval(id_vars(tbl)), .SDcols = col]
    impute_idx <- is.na(tbl[[col]]) & !is.na(tbl[["lag"]])
    tbl[impute_idx, eval(col)] <-  tbl[["lag"]][impute_idx]

  }

  tbl
}

carry_values <- function(tbl, spec) {

  for(s in names(spec)) {
    tbl <- carry_values_single(tbl, s, spec[[s]][["time"]], spec[[s]][["dir"]])
  }

  tbl
}

w_value <- function(source, concept, dir, verbose = FALSE,
  upto = hours(24L), patient_ids = NULL, imp_val = NULL) {


  tbl <- load_concepts(concept, source, patient_ids = patient_ids, verbose = verbose)

  if(is.null(imp_val)) imp_val <- median(tbl[[concept]], na.rm = T)
  if(grepl("epi|ins|norepi", concept)) imp_val <- 0

  if(is_ts_tbl(tbl)) tbl <- tbl[get(index_var(tbl)) <= upto]

  if(is.null(patient_ids)) patient_ids <- unique(stay_windows(source)[[id_vars(tbl)]])

  pts <- data.table::data.table(patient_ids)
  data.table::setnames(pts, names(pts), id_vars(tbl))

  pts <- as_id_tbl(pts)

  w_fun <- ifelse(dir == "increasing", max_or_na, min_or_na)

  res <- tbl[, w_fun(get(concept)), by = eval(id_vars(tbl))]
  data.table::setnames(res, "V1", "w_val")

  res <- merge(pts, res, by = id_vars(res), all.x = T)

  if(length(imp_val) == 2) {
    res[is.na(w_val), "w_val"] <- runif(sum(is.na(res[["w_val"]])), imp_val)
  } else {
    res[is.na(w_val), "w_val"] <- imp_val
  }


  res

}

hypo <- function(source, patient_ids, hypo.threshold = 3.9, upto = hours(72L),
                 verbose = FALSE, value = F) {

  hypo <- load_concepts("glu", source, patient_ids = patient_ids,
                        verbose = verbose)

  hypo[, hg := (glu <= 18.016*hypo.threshold)]
  hypo <- hypo[hg == T & get(index_var(hypo)) <= upto, head(.SD, n = 1L),
               by = eval(id_vars(hypo))]

  if (value) return(hypo[, c(meta_vars(hypo), "glu"), with = FALSE])
  hypo[, c(meta_vars(hypo), "hg"), with = FALSE]

}

fwrap <- function(feat) {
  if(feat == "alt") return("ALT")
  if(feat == "ast") return("AST")
  else return(feat)
}

diff_means <- function(d, ind, threshold) {
  use <- d[ind, ]
  mean(!is.na(use[w_val > threshold, "hg"])) - 
    mean(!is.na(use[w_val <= threshold, "hg"]))
}

ate_sa <- function(d, ind, threshold) {
  use <- d[ind, ]

  H <- as.integer(!is.na(use[["hg"]]))
  L <- as.integer(use[["w_val"]] > threshold)
  S <- use[["w_sofa"]]
  S[is.na(S)] <- 0
  S[S > 12] <- 13

  assert_that(length(unique(S)) == 14L)
  tabular <- table(L, S)
  propensity <- tabular[1, ] / (tabular[2, ] + tabular[1, ])
  mean(H*(2*L-1) / (propensity[S+1]*(1-L) + (1-propensity[S+1])*L))

}

lwrap <- function(lab) {
  if(lab == "liver") return("Liver dysfunction")
  if(lab == "shock_state") return("Shock")
  if(lab == "ins_therapy") return("Insulin therapy")

  lab

}

H_test <- function(x, y, z = NULL) {
  
  if (!is.null(z)) {
    lapply(sort(unique(z)), function(zval) {
      cat("Z-value", zval, "\n")
      H_test(x[z == zval], y[z == zval])
    })
    return()
  }
  
  if (is.logical(y)) {
    cat("Logical, two-sample test with p-value:",
        wilcox.test(x ~ y)$p.value, "\n")
    cat("Logical, contigency test p-value:",
        chisq.test(table(x > 25, y))$p.value, "\n")
  } else{
    cat("Mann-Whitney U-test with p-value:",
        wilcox.test(y[x > 25], y[x <= 25])$p.value, "\n")
  }
  
}

sens_spec_table <- function(score, outcome, src) {
  value_set <- sort(unique(score))
  sens <- spec <- NULL
  for(thresh in value_set) {
    pred <- as.integer(score >= thresh)
    sens <- c(sens, paste0(round(100*sum(pred == 1 & outcome == 1) / 
                                   sum(outcome == 1), 2), "%"))
    spec <- c(spec, paste0(round(100*sum(pred == 0 & outcome == 0) / 
                                   sum(outcome == 0), 2), "%"))
  }

  res <- as.data.frame(cbind(sens, spec), stringsAsFactors = F)
  res <- cbind(threshold = value_set, res)
  names(res) <- c("threshold", paste0(c("sens", "spec"), "_", src))

  res
}

tw_glucose <- function(source, 
                       patient_ids = config("cohort")[[source]][["all"]]) {

  x <- fill_gaps(load_concepts("glu", source, patient_ids = patient_ids, 
                               verbose = F))
  x[, glu := data.table::nafill(glu, "locf"), by = eval(id_var(x))]

  wins <- stay_windows(source)
  x <- merge(x, wins)
  x[get(index_var(x)) >= start & get(index_var(x)) <= end]

  x[, mean(glu, na.rm = T), by = eval(id_var(x))][["V1"]]
}

mean_glucose <- function(source, 
                         patient_ids = config("cohort")[[source]][["all"]]) {

  x <- fill_gaps(load_concepts("glu", source, patient_ids = patient_ids, 
                               verbose = F))
  #x[, glu := data.table::nafill(glu, "locf"), by = eval(id_var(x))]

  wins <- stay_windows(source)
  x <- merge(x, wins)
  x[get(index_var(x)) >= start & get(index_var(x)) <= end]

  x[, mean(glu, na.rm = T), by = eval(id_var(x))][["V1"]]
}

insulin_days <- function(source, 
                         patient_ids = config("cohort")[[source]][["all"]], 
                         upto = hours(10*24)) {

  wins <- stay_windows(source)
  patient_ids <- intersect(id_col(wins[!is.na(end)]), patient_ids)
  wins <- wins[get(id_var(wins)) %in% patient_ids]
  x <- load_concepts("ins", source, patient_ids = patient_ids, verbose = F)


  wins[, end := hours(24*round(end/24))] # round the stay to days
  wins[end > upto, "end"] <- upto
  wins[, num_days := as.integer(end/24)]

  x <- merge(x, wins, all = T)
  x <- x[get(index_var(x)) >= start & get(index_var(x)) < end]

  num_days <- sum(
    x[, length(unique(.bincode(get(index_var(x)), 
                               seq(-0.01, as.integer(upto)-24, 24)))), 
      by = eval(id_var(x))][["V1"]])

  res <- rep(FALSE, sum(wins[["num_days"]]))
  res[1:num_days] <- TRUE

  res
}

remove_doi <- function(src) {
  
  x <- merge(load_concepts("death", src), stay_windows(src), all.x = T)
  id_col(x[get(index_var(x)) < start | get(index_var(x)) > end])
  
}


POC_table <- function(poc, table.names, path) {

  my_doc <- read_docx()

  for (src in names(poc)) {

    res <- Reduce(function(x, y) merge(x, y, by = "group", sort = F),
                  poc[[src]])
    names(res) <- table.names

    my_doc <- body_add(my_doc, srcwrap(src), style = "heading 1")
    my_doc <- my_doc %>% body_add_table(res, style = "table_template")

  }

  print(my_doc, target = path)

  paste("POC table created in", path)

}

PO_char <- function(src, target, breaks = config("bmi-bins")[["who"]], 
                    patient_ids = config("cohort")[[src]][["bmi"]]) {
  
  tbl <- merge(
    get_target(src, target, upto = hours(Inf), patient_ids = patient_ids),
    load_concepts("bmi", src, patient_ids = patient_ids), all.x = TRUE
  )
  
  tbl[, bins := .bincode(bmi, breaks = c(-Inf, breaks, Inf))]
  res <- tbl[, mean(get(target), na.rm = TRUE), by = "bins"]
  res <- data.table::setnames(res, "V1", target)
  res <- data.table::setorderv(res, "bins")
  res[["group"]] <- paste("BMI", bin_labels(breaks, "kg/m2"))
  
  res <- res[, c("group", target), with = F]
  res[[target]] <- round(res[[target]], 2)
  
  res
  
}

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
    
    return(list(val_col, "%", round(100 * sum(x[[val_col]]) / 
                                      length(patient_ids))))
    
  }
  
  list(val_col, "%", round(100 * mean(x[[val_col]])))
  
}

pts_source_sum <- function(source, patient_ids) {
  
  tbl_list <- lapply(
    vars,
    function(x) x[["callback"]](
      load_concepts(x[["concept"]], source, patient_ids = patient_ids, 
                    keep_components = T), unlist(patient_ids)
    )
  )
  
  pts_tbl <- Reduce(rbind,
                    lapply(
                      tbl_list,
                      function(x) data.frame(Reduce(cbind, x))
                    )
  )
  
  cohort_info <- as.data.frame(cbind("Cohort size", "n", 
                                     length(unlist(patient_ids))))
  names(cohort_info) <- names(pts_tbl)
  
  pts_tbl <- rbind(
    cohort_info,
    pts_tbl
  )
  
  names(pts_tbl) <- c("Variable", "Reported", 
                      paste(srcwrap(source), collapse = "-"))
  
  pts_tbl$Variable <- mapvalues(pts_tbl$Variable,
                                from = names(concept_translator),
                                to = sapply(names(concept_translator), 
                                            function(x) concept_translator[[x]])
  )
  
  pts_tbl
  
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

pol_varnames <- function(x) {
  
  subs <- list(
    list("bmi", "BMI"),
    list("glu", "Blood glucose"),
    list("lact", "Blood lactate"),
    list("ins_ifx", "Insulin"),
    list("shock_yes", "MAP < 60 mmHg or vasopressor therapy"),
    list("shock_no", "MAP ≥ 60 mmHg, no vasopressor therapy"),
    list("sofa_cns_comp", "SOFA CNS"),
    list("sofa_coag_comp", "SOFA Coagulation"),
    list("sofa_renal_comp", "SOFA Renal"),
    list("sofa_resp_comp", "SOFA Respiratory"),
    list("DM", "Diabetes"),
    list("source", ""),
    list("cortico", "Corticosteroids"),
    list("enteral", "Enteral nutrition"),
    list("TPN", "Parenteral nutrition"),
    list("dex_amount", "Dextrose 10%"),
    list("sofa_wo_cardio", "SOFA*")
  )
  
  for (i in seq_len(length(subs))) 
    x <- gsub(subs[[i]][[1]], subs[[i]][[2]], x)
  
  adds <- list(
    list("lactate", "mmol/L"),
    list("glucose", "mg/dL"),
    list("BMI", "kg/m2"),
    list("Insulin", "u/h"),
    list("Dextrose", "mL/h")
  )
  
  for (i in seq_len(length(adds))) 
    x <- ifelse(grepl(adds[[i]][[1]], x), paste(x, adds[[i]][[2]]), x)
  
  x
}
