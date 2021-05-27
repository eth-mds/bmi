library(data.table)
library(ggplot2)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

n_orboot <- 100

res <- lapply(
  seq_len(n_orboot), 
  function(i) {
    read.csv(file.path(root, "bootcsv", paste0("boot_", i, ".csv")))       
  }
)

res <- as.data.table(Reduce(rbind, res))

# smoothing
{
  res[, doxl := shift(dox, 1L), by = c("dataset", "job_id")]
  res[, doxr := shift(dox, -1L), by = c("dataset", "job_id")]
  
  res <- cbind(res, doxx = rowMeans(res[, c("dox", "doxl", "doxr")], 
                                    na.rm = TRUE))
  res[bmi == 25L, "doxx"] <- 0
}


plt <- res[, list(Qlow = quantile(doxx, 0.05), 
                  QMid = mean(doxx), 
                  Qhigh = quantile(doxx, 0.95)),
          by = c("dataset", "bmi")]

ggplot() + 
  geom_rect(aes(xmin = 18.5, xmax = 25, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.2, color = NA) +
  geom_ribbon(
    data = plt, 
    aes(x = bmi, ymin = Qlow, ymax = Qhigh, fill = srcwrap(dataset)),
    alpha = 0.2) +
  geom_line(
    data = plt, 
    aes(x = bmi, y = QMid, color = srcwrap(dataset), fill = srcwrap(dataset)),
    size = 2, alpha = 1) + 
  scale_fill_discrete(name = "Dataset") +
  scale_color_discrete(name = "Dataset") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + theme(
    legend.position = c(0.7, 0.7),
    legend.box.background = element_rect()
  ) + xlab("BMI (kg/m2)") + ylab("log(Odds Ratio)") + 
  ggtitle("Non-parametric adjustment of hypoglycemia risk") +
  geom_text(
    aes(x = 22, y = 3, label = "normal BMI range"),
    size = 4#, vjust = 0, hjust = 0, nudge_x = 50, check_overlap = TRUE
  )

ggsave(file.path(root, "figures", "eFigure4.tiff"),
       width = 10, height = 8)
  