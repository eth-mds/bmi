library(ggplot2)
library(magick)
library(ggtext)

feat <- c("Mean glucose [mg/dl]", "Mortality [%]", "Hypoglycemia [%]")

df <- data.frame(
  Group = rep(1:3, 3),
  Value = c(130, 137, 148, 15, 10, 8, 20, 15, 12),
  Feature = factor(rep(feat, each = 3), levels = feat)
)

labs <- c(
  "<b>UNDERWEIGHT</b><br><br><span style='font-size:14pt'>&lt; 18.5</span>",
  "<b>NORMAL</b><br><br><span style='font-size:14pt'>18.5 - 25</span>",
  "<b>OVERWEIGHT</b><br><br><span style='font-size:14pt'>&gt; 25</span>"
)

grid <- ggplot(df, aes(x = Group, y = Value)) +
  geom_line(color = NA) +
  facet_grid(rows = vars(Feature), scales = "free", switch = "y") +
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
  scale_x_continuous(breaks = 1:3, limits = c(0.75, 3.3), labels = labs) +
  scale_y_continuous(position = "right") +
  ggtitle("Hypoglycaemia and the obesity paradox")

plot <- ggplot(df, aes(x = Group, y = Value)) +
  geom_line() +
  facet_grid(rows = vars(Feature), scales = "free", switch = "y") +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.title.y = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.background = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_markdown(size = 14),
    axis.text.x = element_markdown(
      color = "white", fill = c("#91bdde", "#6ac877", "#f78344"),
      margin = unit(c(10, 0, 5, 0), "pt"),
      padding = unit(c(10, 10, 5, 10), "pt")
    )
  ) +
  scale_x_continuous("BMI group [kg/m<sup>2</sup>]", breaks = 1:3,
                     limits = c(0.75, 3.3), labels = labs) +
  scale_y_continuous(position = "right") +
  ggtitle("Hypoglycaemia and the obesity paradox")


ggsave("plot.png", plot, width = 7, height = 5.5, dpi = 300, bg = "transparent")
ggsave("grid.png", grid, width = 7, height = 5.5, dpi = 300, bg = "transparent")

plot <- image_read("plot.png")
grid <- image_read("grid.png")

# export area: x (0, 1100), y (0, 750), size: 2000 X 1364 @ 175 dpi
bg <- image_read("bmi-background.png")
res <- image_mosaic(c(grid, bg, plot))

image_browse(res)
