library(ricu)
library(assertthat)
library(data.table)
library(ranger)
library(magrittr)
library(officer)
library(icd)
library(matrixStats)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config", "dict"))

multivariate <- function(src, coh, diabetes = FALSE) {
  
  dat_name <- paste0(paste0(src, collapse = "_"), coh, diabetes, collapse = "_")
  dat_name <- paste0(dat_name, ".RData")
  vars <- list(glu = list(time = 24L, imp_val = NA_real_),
               lact = list(time = 24L, imp_val = 1), 
               ins_ifx = list(time = 12L, imp_val = 0),
               TPN = list(time = 12L, imp_val = 0),
               enteral = list(time = 12L, imp_val = 0),
               cortico = list(time = 12L, imp_val = 0),
               dex_amount = list(time = 12L, imp_val = 0),
               shock = list(time = 24L, imp_val = 0), 
               bmi = list(time = 0L, imp_val = NA))
  if (diabetes) vars[["DM"]] <- list(time = 0L, imp_val = NA)
  
  # Step 1: get raw data
  if (file.exists(file.path(root, "data", dat_name))) {
    load(file.path(root, "data", dat_name))
  } else {
    pids <- lapply(src, function(x) config("cohort")[[x]][[coh]])
    names(pids) <- src
    if (Reduce(sum, lapply(pids, length)) == 0) pids <- NULL
    glyc <- glycemia_treatment(src, vars = vars, patient_ids = pids, 
                               fill_na = TRUE, sofa = TRUE)
    save(glyc, file = file.path(root, "data", dat_name))
  }
  
  # Step 2: get the design matrix
  rsh <- collect_hypo_cases(glyc)

  # keep relevant cases
  rsh <- rsh[!is.na(glu)] # which proportion is this?!
  rsh <- rsh[glu > 18.016*3.9] # does this throw away anything?!
  reshape_list <- list(
    lact = c(2, 5), glu = 18.016*c(6, 8, 10), ins_ifx = c(0.0001, 2.5, 5),
    shock = 0.5, dex_amount = c(0.0001, 25), 
    bmi = config("bmi-bins")[["who"]]
  )
  
  rsh <- data.frame(rsh)
  rsh <- reshape_frame(rsh, reshape_list)
  
  # Step 3: get the GLM object
  rsh$bmi <- relevel(rsh$bmi, ref = 2)
  if (is.element("source", names(rsh))) {
    rsh$source <- as.factor(rsh$source)
    if (is.element("mimic_demo", rsh$source)) {
      ref_lvl <- "mimic_demo"
    } else if (is.element("aumc", rsh$source)) {
      ref_lvl <- "aumc"
    } else if (is.element("mimic", rsh$source)) {
      ref_lvl <- "mimic"
    }
    rsh$source <- relevel(rsh$source, ref = ref_lvl)
  }
  
  logit0 <- glm(hypo ~ . - bmi, family = "binomial", data = rsh)
  logit <- glm(hypo ~ ., family = "binomial", data = rsh)
  # logit_bmixdiab <- glm(hypo ~ . + bmi:DM, family = "binomial", data = rsh)

  print(anova(logit0, logit, test = "Chisq")$`Pr(>Chi)`[2])
  # print(anova(logit, logit_bmixdiab, test = "Chisq")$`Pr(>Chi)`[2])
 
  # Step 4: get the table with ORs and CIs
  coef <- summary(logit)$coefficients
  
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
    names(logit$xlevels), sapply(logit$xlevels, function(fct) fct[1L])
  )
  
  res <- rbind(
    res,
    cbind(
      ref_names, "1\\n (-)"
    ), use.names = F
  )
  res$V1 <- pol_varnames(res$V1)
  
  res
  
}

all_src <- c("aumc", "hirid", "mimic", "eicu")
sens_src <- c("aumc", "hirid")
us_src <- c("mimic", "eicu")

tbl <- Reduce(
  function(x, y) merge(x, y, by = "V1", all = TRUE),
  list(
    multivariate(all_src, "insulin", diabetes = FALSE),
    multivariate(sens_src, "new", diabetes = FALSE),
    multivariate(us_src, "insulin", diabetes = TRUE)
  )
)

names(tbl) <- c("Odds Ratio (95% CI)", "All datasets", "Sensitivity Analysis", 
                "Diabetes adjusted (US datasets)")


if (all(c("aumc", "hirid", "mimic", "eicu") %in% all_src)) {
  
  my_doc <- read_docx()
  
  my_doc <- my_doc %>%
    body_add_table(tbl, style = "table_template")
  
  print(my_doc, target = file.path(root, "tables", "Table4.docx"))
  
}