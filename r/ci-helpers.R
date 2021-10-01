
CI_dat <- function(src, upto = hours(Inf), y = "death", x = "bmi", 
                   x_bins = config("bmi-bins")[["who"]], 
                   z = NULL, z_binning = NULL, coh = "bmi",
                   subset_fn = NULL,
                   patient_ids = config("cohort")[[src]][[coh]],
                   add_prop = NULL) {
  
  cat("Dataset(s):", src, "\n")
  if (length(src) > 1L) {
    
    patient_ids <- lapply(src, function(dsrc) config("cohort")[[dsrc]][[coh]])
    names(patient_ids) <- src
    
  }
  if (!is.null(subset_fn)) patient_ids <- subset_fn(src, patient_ids, upto)
  
  yt <- get_target(src, y, upto, patient_ids)
  
  # load x
  by.args <- x
  xt <- load_concepts(x, src, patient_ids = patient_ids, verbose = FALSE)
  
  # apply x_bins to xt
  xt[, raw := get(x)]
  xt[[x]] <- .bincode(xt[[x]], c(-Inf, x_bins, Inf))
  
  if (y == "fhm") {
    
    yt <- merge(yt, xt, all.x = TRUE)
    
  } else yt <- merge(yt, xt, all = TRUE)
  
  if (is.logical(yt[[y]])) yt[is.na(get(y)), y] <- FALSE
  
  # load and apply z_bins to z if exists
  if (!is.null(z)) {

    zt <- load_concepts(z, src, patient_ids = patient_ids)
    if (is_ts_tbl(zt)) {
      zt <- zt[, mean(get(z)), by = c(id_var(zt))]
      zt <- rename_cols(zt, z, "V1")
    }
    by.args <- c(by.args, z)
    
    zt[[z]] <- z_binning(zt[[z]])
    
    yt <- merge(yt, zt, all.x = TRUE)
    if (is.logical(zt[[z]])) yt[is.na(get(z)), z] <- FALSE
    
    yt <- yt[!is.na(get(z))]
  }

  # independence tests:
  H_test(x = yt[["raw"]], y = yt[[y]], z = if (!is.null(z)) yt[[z]] else NULL) 
  
  # bootstrap the mean of y and ci(y)
  get_ci <- function(sample) {
    
    if(all(sample == mean(sample))) return(c(mean(sample), mean(sample)))
    sample <- sample[!is.na(sample)]
    boot.samp <- boot(data = sample, statistic =
                        function(data, indices) mean(data[indices], na.rm = F),
                      R = 500)
    boot.ci <- boot.ci(boot.samp, type = "basic")
    
    return(boot.ci$basic[4:5])
    
  }

  ret <- yt[, list(mean(get(y), na.rm = T), list(get_ci(get(y))),
                   Feature = x), by = by.args]
  
  ret[["lower"]] <- sapply(ret[["V2"]], `[[`, 1L)
  ret[["upper"]] <- sapply(ret[["V2"]], `[[`, 2L)
  if (any(ret[["lower"]] < 0)) {
    cat("Concept", y, "has a CI with lower end < 0; Fixing this\n")
    ret[lower < 0, lower := 0]
  }
  
  ret <- ret[, c(x, "V1", "lower", "upper", "Feature", z), with=FALSE]
  data.table::setnames(ret, x, "meanval")
  
  if (!is.null(add_prop)) {
    pat_prop <- 100 * length(unique(id_col(yt[get(y) >0]))) / 
      length(patient_ids)
    ret[, dataset := paste0(srcwrap(src), " (", spec_dec(pat_prop, 1),"%)")]
  } else ret[, dataset := paste(src, collapse = ", ")]
  
  ret
  
}

CI_plot <- function(ret, x_label = "BMI (kg/m^2)", y_label = "Mortality", 
                    title = "Mortality by BMI group", z_cond = NULL,
                    z_cond_name = NULL, z_cond_labels = NULL,
                    x_bins = config("bmi-bins")[["who"]],
                    pos.x = 0.7, pos.y = 0.7) {
  
  if (is.data.table(ret) & !is.null(z_cond)) {
    
    ret <- data.table::setnames(ret, z_cond, "z_cond")
    
    p <- ggplot(ret, aes(x = meanval, y = V1, color = factor(z_cond))) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill = factor(z_cond)), 
                  alpha=0.2, size = 0) +
      theme_bw(15) + xlab(x_label) + ylab(y_label) +
      ggtitle(title) + 
      scale_x_continuous(labels=bin_labels(x_bins, NULL),
                         breaks = c(1:(length(x_bins)+1))) + 
      scale_color_discrete(name = z_cond_name, labels = z_cond_labels) +
      scale_fill_discrete(name = z_cond_name, labels = z_cond_labels) +
      theme(
        legend.box.background = element_rect(),
        legend.position = c(pos.x, pos.y)
      )
    
    return(p)
    
  } else if (is.data.table(ret)) {
    
    p <- ggplot(ret, aes(x = meanval, y = V1, color = Feature)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin=lower, ymax=upper, color = Feature), alpha=0.2) +
      theme_bw(15) + xlab(x_label) + ylab(y_label) +
      ggtitle(paste0(srcwrap(src))) +
      scale_x_continuous(labels=bin_labels(x_bins, NULL),
                         breaks = c(1:(length(x_bins)+1)))
    
    return(p)
    
  } else {
    
    ret <- Reduce(rbind, ret)
    
    p <- ggplot(ret, aes(x = meanval, y = V1, color = srcwrap(dataset))) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill = srcwrap(dataset)), 
                  alpha=0.2, size = 0) +
      theme_bw(15) + xlab(x_label) + ylab(y_label) +
      ggtitle(title) +
      scale_x_continuous(labels=bin_labels(x_bins, NULL),
                         breaks = c(1:(length(x_bins)+1))) +
      scale_color_discrete(name = "Dataset") + 
      scale_fill_discrete(name = "Dataset") +
      theme(
        legend.position = c(pos.x, pos.y),
        legend.box.background = element_rect()
      )
    
    return(p)
    
  }
  
}
