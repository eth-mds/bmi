library(ggplot2)
library(reshape2)

miss_pattern <- function(src, tbl) {

  tbl <- merge(tbl, load_concepts("version", src), all.x = TRUE)
  
  vars <- c("glu", "lact", "shock", "bmi")
  
  # patient level missingness
  get_plm <- function(tbl, vars) {
    plm <- tbl[, lapply(.SD, function(x) all(is.na(x))), 
               by = c(id_vars(tbl), "version"),
               .SDcols = vars][, c(vars, "version"), with = FALSE]
    cbind(plm[, lapply(.SD, mean), by = "version"], dataset = src, 
          type = "Patient")
  }
  
  # time-point level missingness
  get_tlm <- function(tbl, vars) {
    tbl <- tbl[get(index_var(tbl)) >= hours(6L) & 
                 get(index_var(tbl)) <= hours(240L)]
    tbl <- tbl[as.integer(get(index_var(tbl))) %% 6 == 0]

    tlm <- tbl[, lapply(.SD, function(x) all(is.na(x))), 
               by = c(meta_vars(tbl), "version"),
               .SDcols = vars][, c(vars, "version"), with = FALSE]
    
    cbind(tlm[, lapply(.SD, mean), by = "version"], dataset = src, 
          type = "Time-point")
  }
  
  rbind(get_plm(tbl, vars), get_tlm(tbl, vars))
  
}

src <- c("mimic", "hirid", "aumc")
dat <- lapply(
  src,
  function(dsrc) {
    glycemia_treatment(dsrc, patient_ids = id_col(
      load_concepts("age", dsrc)[age >= 18L]))
  }
)

mp <- Reduce(rbind, Map(miss_pattern, src, dat))

ggplot(melt(mp), aes(x = variable, y = value, fill = type)) +
  geom_col(position = "dodge") + theme_bw() + ylim(c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(rows = vars(srcwrap(dataset), version))
