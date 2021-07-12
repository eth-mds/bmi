library(ricu)

eicu_search <- function(names, table, col_name, ...) {
  strn <- eicu[[table]][[col_name]]
  a <- lapply(names, function(nm) {
    rows <- grep(nm, strn, ...)
    if(length(rows) == 0) return(NULL)
    res <- table(strn[rows])
    ret <- data.frame(
      item_name = attr(res, "dimnames")[[1]],
      count = as.integer(res)
    )
    ret <- ret[order(-ret$count),]
    return(ret)
  })
  names(a) <- names
  return(a)
}

mimic_search <- function(names, ...) {
  ret <- lapply(names, function(p) {
    rows <- grep(p, mimic[["d_items"]][["label"]], ...)
    if(length(rows) == 0) return(NULL)
    res <- lapply(rows, function(i) {
      item <- mimic[["d_items"]][["itemid"]][i]
      name <- mimic[["d_items"]][["label"]][i]
      table <- mimic[["d_items"]][["linksto"]][i]
      if (table == "microbiologyevents"){
        count <- 0
      } else count <- nrow(subset(mimic[[table]], itemid == item))
      return(c(item, name, table, count))
    })
    res <- Reduce(rbind, res)
    res <- data.frame(matrix(res, ncol = 4))
    if(length(names(res)) != 4) browser()
    names(res) <- c("item", "name", "table", "count")
    res$item <- as.integer(as.character(res$item))
    res$count <- as.integer(as.character(res$count))
    res <- res[order(-res$count),]
    rownames(res) <- NULL
    return(res)
  })
  names(ret) <- names
  return(ret)
}

miiv_search <- function(names) {
  ret <- lapply(names, function(p) {
    rows <- grep(p, miiv[["d_items"]][["label"]], ignore.case = TRUE)
    if(length(rows) == 0) return(NULL)
    res <- lapply(rows, function(i) {
      item <- miiv[["d_items"]][["itemid"]][i]
      name <- miiv[["d_items"]][["label"]][i]
      table <- miiv[["d_items"]][["linksto"]][i]
      if (table == "microbiologyevents"){
        count <- 0
      } else count <- nrow(subset(miiv[[table]], itemid == item))
      return(c(item, name, table, count))
    })
    res <- Reduce(rbind, res)
    res <- data.frame(matrix(res, ncol = 4))
    if(length(names(res)) != 4) browser()
    names(res) <- c("item", "name", "table", "count")
    res$item <- as.integer(as.character(res$item))
    res$count <- as.integer(as.character(res$count))
    res <- res[order(-res$count),]
    rownames(res) <- NULL
    return(res)
  })
  names(ret) <- names
  return(ret)
}

aumc_search <- function(names, tbl, ...) {
  
  ret <- lapply(
    names, function(nm) {
      idx <- grepl(nm, aumc[[tbl]]$item, ...)
      res <- data.table::as.data.table(aumc[[tbl]][idx, ])
      res <- res[, c("itemid", "item"), with = FALSE]
      res <- res[, .N, by = c("itemid", "item")]
      data.table::setorderv(res, "N", -1L)
      
      unique(res)
    }
  )
  names(ret) <- names

  ret
}

us_cort <- c("hydrocortison", "Cortef", "cortison", "ethamethasoneb", 
             "Celeston", "prednison", "intensol", "prednisolon", "Orapred", 
             "Prelone", "triamcinolon", "Kenalog", "Aristospan", "Medrol", 
             "dexamethason", "dexame", "DexPak", "florinef", "cort")

eicu_search(us_cort, "infusiondrug", "drugname")
eicu_search(us_cort, "medication", "drugname")
mimic_search(us_cort)
miiv_search(us_cort)
aumc_search(us_cort, "drugitems", ignore.case = TRUE)

eicu_regex <- c("HYDROCORTISONE", "PREDNISONE", "METHYLPREDNISOLONE",
                "solumedrol", "solu-medrol", "DEXAMETHASONE")

lapply(
  us_cort, 
  function(crt) {
    sort(table(grep(crt, mimic$prescriptions$drug, 
                    ignore.case = TRUE, value = TRUE)),
         decreasing = TRUE)
  }
)

tpn <- mimic_search(c("PN"))
mimic_search("enteral", ignore.case = TRUE)
mimic_search("nutri", ignore.case = TRUE)

# 3428, 3429, 3430, 3427, 5908 -> feed tube!
# inputevents_cv: 30090
load_concepts(c("enteral", "dbsource"), "mimic")
res <- merge(load_concepts("enteral", "mimic"), 
             load_concepts("dbsource", "mimic"))
table(res$dbsource)

eicu_search("PN", "infusiondrug", "drugname")
eicu_search("PN", "medication", "drugname")
eicu_search("PN", "intakeoutput", "celllabel")

###
tpn <- c(
  "trophamin", "tralement", "SMOFlipid", "renamin", "prosol", "procalamin",
  "premasol", "plenamine", "plasma-lyte", "perikabiven", "peditrace", "omegaven",
  "nutrilyte", "nutrilipid", "novamin", "kabiven", "hepatamin", "freamin", 
  "branchamin", "aminosyn", "aminoprotect", "addamel", "freka", "smofkabiven",
  "structokabiven", "nutriflex", "tpn", "parenteral"
)

mimic_search(tpn, ignore.case = TRUE)
lapply(
  tpn, 
  function(crt) {
    sort(table(grep(crt, mimic$prescriptions$drug, 
                    ignore.case = TRUE, value = TRUE)),
         decreasing = TRUE)
  }
)

enteral <- c(
  "novasource", "promote", "isosource", "survimed", "nepro",
  "feed", "nutri", "enteral", "tube"
)

mimic_search(enteral, ignore.case = TRUE)
mimic_search("tube", ignore.case = TRUE)
lapply(
  enteral, 
  function(crt) {
    sort(table(grep(crt, mimic$prescriptions$drug, 
                    ignore.case = TRUE, value = TRUE)),
         decreasing = TRUE)
  }
)

eicu_search(tpn, "infusiondrug", "drugname", ignore.case = TRUE)
eicu_search(tpn, "medication", "drugname", ignore.case = TRUE)
eicu_search(tpn, "intakeoutput", "celllabel", ignore.case = TRUE)

eicu_search(enteral, "infusiondrug", "drugname", ignore.case = TRUE)
eicu_search(enteral, "medication", "drugname", ignore.case = TRUE)
eicu_search(enteral, "intakeoutput", "celllabel", ignore.case = TRUE)

### AUMC ###
aumc_search("tube", "processitems", ignore.case = TRUE)
aumc_search("tube", "procedureorderitems", ignore.case = TRUE)

aumc_search(enteral, "drugitems", ignore.case = TRUE)
aumc_search(tpn, "drugitems", ignore.case = TRUE)
aumc_search("voeding", "drugitems", ignore.case = TRUE)

aumc_med <- as.data.table(aumc$drugitems)
aumc_med <- aumc_med[, .N, by = c("item", "itemid")]
setorderv(aumc_med, "itemid")

aumc_med[itemid %in% (-100:100 + 20730)]

library(ricu)
tbl <- load_concepts("enteral", "mimic")
quantile(tbl$dur_var, seq(0.1, 1, 0.1))

load_concepts("enteral", "hirid")
load_concepts("enteral", "eicu")


load_concepts("TPN", src)
load_concepts("cortico", "mimic", 
              patient_ids = config("cohort")[["mimic"]][["insulin"]])



