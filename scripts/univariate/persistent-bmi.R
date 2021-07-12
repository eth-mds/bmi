library(ricu)
library(ggplot2)
library(assertthat)
library(data.table)

invisible(lapply(list.files(here::here("r"), full.names = TRUE), 
                 source))

persistent_bmi <- function(src = c("mimic", "eicu", "hirid", "aumc"),
                           los_thresh = c(0, 4, 6, 8, 10),
                           breaks = c(-Inf, 
                                      config("bmi-bins")[["who"]], 
                                      Inf)
                           ) {
  
  tbl <- load_concepts(c("bmi", "los_icu", "death"), src)
  tbl <- tbl[!is.na(bmi)]
  tbl[, bmi_bin := .bincode(bmi, breaks)]
  
  res <- lapply(
    los_thresh,
    function(th) {
      ret <- tbl[los_icu > th]
      ret <- ret[, list(death = any(death), bmi_bin = bmi_bin), 
          by = c(id_vars(ret))]
      ret[is.na(death), death := FALSE] 
      ret <- ret[, list(death_rate = mean(death)), by = "bmi_bin"]
      ret[, los_cutoff := th]
    }
  )
  
  Reduce(rbind, res)
  
}

res <- persistent_bmi()
ggplot(res, aes(x = bmi_bin, y = death_rate, 
                color = factor(los_cutoff))) +
  geom_line() + theme_bw() +
  xlab("BMI group (kg/m^2)") + ylab("Mortality rate") +
  scale_x_continuous(
    labels = bin_labels(config("bmi-bins")[["who"]], NULL),
    breaks = seq.int(1L, 6L, 1L)
  ) +
  scale_color_discrete(name = "ICU LOS > than (days)") +
  theme(
    legend.position = c(0.7, 0.7),
    legend.box.background = element_rect()
  )
