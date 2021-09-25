
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

library(ggplot2)
library(magick)
library(ggtext)

feat <- c(tw_avg_glucose = "Mean glucose [mg/dl]", death = "Mortality [%]",
          hypo = "Hypoglycemia [%]")
bins <- c(18.5, 25, 30)

res <- Map(CI_dat, y = names(feat), MoreArgs = list(
  src = c("aumc", "hirid", "mimic", "eicu"),
  x_bins = bins,
  patient_ids = lapply(config("cohort"), `[[`, "bmi")
))

res <- Map(`cbind`, res, Name = lapply(factor(feat, levels = feat), rep,
                                       length(bins) + 1L))
res <- rbind_lst(res)
res <- rename_cols(res[, list(meanval, V1, Name)],
                   c("Group", "Value", "Feature"))

res$Value <- ifelse(grepl("\\[%\\]$", res$Feature), res$Value * 100, res$Value)

res <- cbind(as.data.frame(res),
             Ymin = rep(c(131.5, 7, 11.5), each = 4),
             Ymax = rep(c(150, 14, 23), each = 4))

labs <- paste0(
  "<b>",
  c("UNDERWEIGHT", "NORMAL", "OVERWEIGHT", "OBESE (CLASS I-III)"),
  "</b><br><br><span style='font-size:14pt'>",
  c(paste("&lt;", bins[1L]), paste(bins[1L], "-", bins[2L]),
    paste(bins[2L], "-", bins[3L]), paste("&gt;", bins[3L])),
  "</span>"
)

grid <- ggplot(res, aes(x = Group, y = Value)) +
  geom_line(color = NA) +
  facet_grid(rows = vars(Feature), scales = "free_y", switch = "y") +
  geom_blank(aes(y = Ymin)) +
  geom_blank(aes(y = Ymax)) +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.title = element_blank(), axis.ticks = element_blank(),
    axis.text = element_text(color = NA), strip.text = element_text(color = NA),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.title.x = element_markdown(size = 14, color = NA),
    axis.text.x = element_markdown(
      color = NA, fill = NA,
      margin = unit(c(10, 0, 5, 0), "pt"),
      padding = unit(c(10, 10, 5, 10), "pt")
    ), plot.title = element_text(color = NA)
  ) +
  scale_x_continuous(breaks = 1:4, limits = c(0.75, 4.35), labels = labs) +
  scale_y_continuous(position = "right") +
  ggtitle("Hypoglycemia and the obesity paradox")

plot <- ggplot(res, aes(x = Group, y = Value)) +
  geom_line() +
  facet_grid(rows = vars(Feature), scales = "free_y", switch = "y") +
  geom_blank(aes(y = Ymin)) +
  geom_blank(aes(y = Ymax)) +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.title.y = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.background = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_markdown(size = 14),
    axis.text.x = element_markdown(
      color = "white", fill = c("#91bdde", "#6ac877", "#e3e23e", "#f78344"),
      margin = unit(c(10, 0, 5, 0), "pt"),
      padding = unit(c(10, 10, 5, 10), "pt")
    )
  ) +
  scale_x_continuous("BMI group [kg/m<sup>2</sup>]", breaks = 1:4,
                     limits = c(0.75, 4.35), labels = labs) +
  scale_y_continuous(position = "right") +
  ggtitle("Hypoglycemia and the obesity paradox")

ggsave(file.path(root, "scripts", "graph-abstract", "plot.png"), plot,
       width = 8, height = 5.5, dpi = 300, bg = "transparent")
ggsave(file.path(root, "scripts", "graph-abstract", "grid.png"), grid,
       width = 8, height = 5.5, dpi = 300, bg = "transparent")

plot <- image_read(file.path(root, "scripts", "graph-abstract", "plot.png"))
grid <- image_read(file.path(root, "scripts", "graph-abstract", "grid.png"))

# export area: x (0, 1100), y (0, 750), size: 2000 X 1364 @ 175 dpi
bg <- image_read(file.path(root, "scripts", "graph-abstract",
                           "bmi-background.png"))
res <- image_mosaic(c(grid, bg, plot))

image_browse(res)

image_write(image_resize(res, "1200x"),
            path = file.path(root, "scripts", "graph-abstract", "graph-abstract.png"),
            format = "png")
