library(ricu)
library(ggplot2)
library(assertthat)
library(boot)
library(data.table)
library(scales)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

Sys.setenv(RICU_CONFIG_PATH = file.path(root, "config", "dict"))

src <- c("aumc", "hirid", "mimic", "eicu") # c("mimic_demo", "eicu_demo") # 
cnc <- c("crea", "lact", "bili", "plt", "pafi", "map", "gcs", "glu")
coh <- lapply(src, function(x) config("cohort")[[x]][["insulin"]])
names(coh) <- src
lwr <- hours(0L)
upr <- hours(24L)

data <- list()
for (dsrc in src) {
  
  tbl <- load_concepts(c(cnc, "bmi"), dsrc, 
                       patient_ids = config("cohort")[[dsrc]][["bmi"]])
  
  tbl[, source := dsrc]
  tbl <- tbl[get(index_var(tbl)) >= lwr & get(index_var(tbl)) <= upr]
  #browser()
  tbl <- as_ts_tbl(tbl, id_vars = c(id_var(tbl), "source"))
  
  data[[dsrc]] <- tbl
  
}

tbl <- Reduce(rbind, data)
tbl[, bmi := .bincode(bmi, c(-Inf, config("bmi-bins")[["who"]], Inf))]

tbl <- tbl[, lapply(.SD, function(x) all(is.na(x))), by = c(id_vars(tbl), "bmi"), 
           .SDcols = cnc]

res <- tbl[, lapply(.SD, mean), by = c("source", "bmi"), 
    .SDcols = cnc]

setorderv(res, c("source", "bmi"))

res[, source := srcwrap(source)]
setnames(res, "source", "Dataset")
plt <- list()
fig <- list(
  crea = list(full = "Creatinine", pos = c(0.7, 0.3)),
  lact = list(full = "Lactate", pos = c(0.8, 0.7)),
  bili = list(full = "Bilirubin", pos = c(0.6, 0.2)),
  plt = list(full = "Platelets", pos = c(0.3, 0.4)),
  pafi = list(full = "PaO2/FiO2", pos = c(0.7, 0.7)),
  map = list(full = "MAP", pos = c(0.7, 0.8)),
  gcs = list(full = "Glasgow Coma Scale", pos = c(0.6, 0.4)),
  glu = list(full = "Glucose", pos = c(0.3, 0.5))
)

for (c in cnc) {
  
  plt[[c]] <- ggplot(res, aes_string(x = "bmi", y = c, color = "Dataset")) +
    geom_line() + theme_bw() +
    ggtitle(paste0(fig[[c]][["full"]], " reporting at 24 hours")) +
    theme(
      legend.position = fig[[c]][["pos"]],
      legend.box.background = element_rect()
    ) + scale_y_continuous(labels = percent) +
    #scale_x_continuous(labels=bin_labels(config("bmi-bins")[["who"]], NULL)) +
    ylab("Percent patients with missing value") + xlab("BMI (kg/m2)")
  
}

# plot
cowplot::plot_grid(plotlist = plt, ncol = 4L)
ggsave(file.path(root, "4files-BMI", "eSupplement", "eFigure3.tiff"),
       width = 20, height = 10)

# Table
mkt <- data.table::copy(res)
mkt[, bmi := bin_labels(config("bmi-bins")[["who"]], NULL)[bmi]]
mkt <- mkt[, lapply(.SD, function(x) round(100 * x)), by = c("bmi", "Dataset")]
setnames(mkt, cnc, vapply(fig, function(x) x[["full"]], character(1L)))
setnames(mkt, "bmi", "BMI group")
setcolorder(mkt, "Dataset")

my_doc <- read_docx()
my_doc <- my_doc %>% body_add_table(mkt, style = "table_template")
print(my_doc, target = file.path(root, "tables", "eSupplement_Table2.docx"))
