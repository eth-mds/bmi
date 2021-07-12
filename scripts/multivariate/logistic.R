library(ricu)
library(assertthat)
library(data.table)
library(ranger)
library(magrittr)
library(officer)
library(matrixStats)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

load_bsgi <- function(data_source) {
  
  vars <- list(glu = list(time = 24L, imp_val = NA_real_),
               lact = list(time = 24L, imp_val = 1), 
               ins = list(time = 12L, imp_val = 0),
               shock = list(time = 24L, imp_val = 0), 
               bmi = list(time = 0L, imp_val = NA),
               weight = list(time = 0L, imp_val = NA), 
               liver_damage = list(time = 48L, imp_val = 0))
  
  glyc <- glycemia_treatment(data_source, vars = vars,
                             patient_ids = 
                               config("cohort")[[undemo(data_source)]][["insulin"]], 
                             fill_na = T)
  
  glyc[, bmi := nafill(bmi, "locf"), by = eval(id_var(glyc))]
  glyc[, weight := nafill(weight, "locf"), by = eval(id_var(glyc))]
  
  glyc
  
}

model_summary <- function(glyc, ins_norm = F, reshape = TRUE) {
  
  glyc <- copy(glyc)
  if (ins_norm) glyc[, ins := (ins / weight)]
  glyc[, weight := NULL]

  rsh <- collect_hypo_cases(glyc)

  # keep relevant cases
  rsh <- rsh[!is.na(glu)]
  rsh <- rsh[glu > 18.016*3.9]
  reshape_list <- list(
    lact = c(2, 5),
    glu = 18.016*c(6, 8, 10),
    ins = c(0.001, 2.5, 5) / (1 + 74*ins_norm), liver_damage = 0.5,
    shock = 0.5, bmi = config("bmi-bins")[["who"]]
  )
  
  if (reshape) {
    rsh <- data.frame(rsh)
    rsh <- reshape_frame(rsh, reshape_list)
  } 

  rsh

}

src <- c("mimic", "eicu", "hirid", "aumc") # c("mimic_demo", "eicu_demo") #

glyc <- lapply(src, load_bsgi)
res <- lapply(glyc, model_summary)

ref_bmi_normal <- TRUE

res_logit <- lapply(res, function(x) {

  if (ref_bmi_normal) {

    x$bmi <- relevel(x$bmi, ref = 2)

  }

  logit0 <- glm(hypo ~ . - bmi, family = "binomial", data = x)
  logit <- glm(hypo ~ ., family = "binomial", data = x)

  print(anova(logit0, logit, test = "Chisq")$`Pr(>Chi)`[2])

  logit

  }
)

est_orc <- function(x) {

  coef <- summary(x)$coefficients

  OR <- cbind(
    coef[, "Estimate"],
    coef[, "Estimate"] - 1.96 * coef[, "Std. Error"],
    coef[, "Estimate"] + 1.96 * coef[, "Std. Error"]
  )

  OR <- round(exp(OR), 2)


  res <- data.table(cbind(
    rownames(coef),
    paste0(OR[, 1], "\\n (", OR[, 2], "-", OR[, 3], ")")
  ))
  
  # add the baseline bands
  ref_names <- paste0(
    names(x$xlevels), sapply(x$xlevels, function(fct) fct[1L])
  )
  
  #browser()
  
  rbind(
    res,
    cbind(
      ref_names, "1\\n (-)"
    ), use.names = F
  )

}

tbl <- Reduce(
  function(x, y) merge(x, y, by = "V1"),
  lapply(res_logit, est_orc)
)

names(tbl) <- c("Odds Ratio (95% CI)", sapply(src, srcwrap))

tbl[[1]] <- gsub("bmi", "BMI", tbl[[1]])
tbl[[1]] <- gsub("glu", "Blood glucose", tbl[[1]])
tbl[[1]] <- gsub("lact", "Blood lactate", tbl[[1]])
tbl[[1]] <- gsub("ins", "Insulin", tbl[[1]])
tbl[[1]] <- gsub("shock_yes", "MAP < 60 mmHg or vasopressor therapy", tbl[[1]])
tbl[[1]] <- gsub("shock_no", "MAP ≥ 60 mmHg, no vasopressor therapy", tbl[[1]])
tbl[[1]] <- gsub("liver_damage_yes", "ALT/AST ≥ 40 IU or bilirubin ≥ 2.0 mg/dL", 
                 tbl[[1]])
tbl[[1]] <- gsub("liver_damage_no", "ALT/AST < 40 IU, bilirubin < 2.0 mg/dL", 
                 tbl[[1]])

tbl[[1]] <- gsub("sofa_cns_comp", "SOFA CNS", tbl[[1]])
tbl[[1]] <- gsub("sofa_coag_comp", "SOFA Coagulation", tbl[[1]])
tbl[[1]] <- gsub("sofa_renal_comp", "SOFA Renal", tbl[[1]])
tbl[[1]] <- gsub("sofa_resp_comp", "SOFA Respiratory", tbl[[1]])

tbl[[1]] <- ifelse(grepl("lactate", tbl[[1]], ignore.case = T), 
                   paste(tbl[[1]], "mmol/L"), tbl[[1]])
tbl[[1]] <- ifelse(grepl("glucose", tbl[[1]], ignore.case = T), 
                   paste(tbl[[1]], "mg/dL"), tbl[[1]])
tbl[[1]] <- ifelse(grepl("bmi", tbl[[1]], ignore.case = T), 
                   paste(tbl[[1]], "kg/m2"), tbl[[1]])
tbl[[1]] <- ifelse(grepl("insulin", tbl[[1]], ignore.case = T), 
                   paste(tbl[[1]], "u/h"), tbl[[1]])


if (all(c("mimic", "eicu", "hirid", "aumc") %in% src)) {
  
  my_doc <- read_docx()
  
  my_doc <- my_doc %>%
    body_add_table(tbl, style = "table_template")
  
  print(my_doc, target = file.path(root, "tables", "Table2.docx"))
  
}