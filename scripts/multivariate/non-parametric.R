#!/usr/bin/env Rscript

#BSUB -W 00:30
#BSUB -n 12
#BSUB -R rusage[mem=16000]
#BSUB -J boot[1-5]
#BSUB -o bootcsv/boot_%J.out
library(matrixStats)
library(data.table)
library(ranger)
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

jid <- as.integer(Sys.getenv("LSB_JOBINDEX"))
load("bmi/orboot.RData")

cts_adjust <- function(data, src, seed, bmi_seq = seq.int(15, 35, 5)) {
  
  data <- copy(data)
  norm_ind <- which(bmi_seq == 25L)
  folds <- sample.int(nrow(data), round(nrow(data) * 0.8), replace = FALSE)
  
  c_eff <- function(form, data, folds, bmi_seq = seq.int(15, 35, 5), 
                    ...) {
    
    rf <- ranger(form, data = data[folds], keep.inbag = T, 
                 importance = "impurity", 
                 probability = TRUE, seed = seed, num.threads = 12, 
                 ...)
    oob.matrix <- Reduce(cbind, lapply(rf$inbag.counts, function(x=i) x == 0))
    
    p2 <- lapply(
      bmi_seq,
      function(do_bmi) {
        
        data[, bmi := do_bmi]
        prob <- predict(rf, data[folds], predict.all = T,
                        num.threads = 12)$predictions[, 2, ]
        #browser()
        prob <- rowSums(prob * oob.matrix) / rowSums(oob.matrix)
        prob[prob == 0] <- min(prob[prob > 0])
        prob[prob == 1] <- max(prob[prob < 1])
        prob / (1 - prob)
        
      }
    )
    
    p2 <- Reduce(cbind, p2)
    p2 <- p2 / p2[, norm_ind]
    
    colMedians(log(p2))
    
  }
  
  data.frame(
    bmi = bmi_seq,
    dox = c_eff(form = hypo ~ ., data = data, folds = folds,
                bmi_seq = bmi_seq),
    dataset = src,
    job_id = seed
  )
  
}

src <- c("mimic", "eicu", "hirid", "aumc")
or <- Map(cts_adjust, res_cts, src, jid)
or <- Reduce(rbind, or)

write.csv(or, file = paste0("bootcsv/boot_", jid, ".csv"))