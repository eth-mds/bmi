library(ricu)
library(stringr)
library(magrittr)
library(officer)
library(assertthat)
library(plyr)

root <- rprojroot::find_root(".gitignore")
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

src <- c("aumc", "hirid", "mimic", "eicu")

po_chars <- list(
  glu = list(target = "tw_avg_glucose", coh = "bmi"),
  glu_freq = list(target = "glu_freq", coh = "bmi"),
  ins = list(target = "max_insulin", coh = "insulin"),
  ins_wnorm = list(target = "max_insulin_wnorm", coh = "insulin"),
  dur_TPN = list(target = "dur_TPN", coh = "insulin"),
  dur_enteral = list(target = "dur_enteral", coh = "insulin"),
  dur_cortico = list(target = "dur_cortico", coh = "insulin"),
  tw_avg_dextrose = list(target = "tw_avg_dextrose", coh = "insulin")
)

poc <- lapply(
  src,
  function(src) {
    lapply(
      po_chars,
      function(x) PO_char(src, x[["target"]], 
                          patient_ids = config("cohort")[[src]][[x[["coh"]]]])
    )
  }
)
names(poc) <- src

POC_table(poc, c("BMI group", names(po_chars)), 
          file.path(root, "tables", "Table3.docx"))



