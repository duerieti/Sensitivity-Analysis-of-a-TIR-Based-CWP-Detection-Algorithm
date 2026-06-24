rnorm(100,0,1.5)
rnorm(100,4,0.5)
rnorm(100,4,1.5)
rnorm(100,0,0.5)


library(ggplot2)


# Generate data for all four panels
data <- rbind(
  data.frame(
    value = rnorm(200, 0, 1.7),
    panel = "Small Effect with nonlinearity/Interactions"
  ),
  data.frame(
    value = rnorm(200, 4, 0.7),
    panel = "Large Effect with no nonlinearity/Interactions"
  ),
  data.frame(
    value = rnorm(200, 4, 1.7),
    panel = "Large Effect with nonlinearity/Interactions"
  ),
  data.frame(
    value = rnorm(200, 0, 0.7),
    panel = "Small Effect with no nonlinearity/Interactions"
  )
)

data$panel <- factor(data$panel, levels = c(
  "Small Effect with nonlinearity/Interactions",
  "Large Effect with nonlinearity/Interactions",
  "Small Effect with no nonlinearity/Interactions",
  "Large Effect with no nonlinearity/Interactions"
))

# Compute per-panel statistics for annotations
panel_stats <- do.call(rbind, lapply(levels(data$panel), function(p) {
  x <- data$value[data$panel == p]
  data.frame(panel = factor(p, levels = levels(data$panel)),
             mu = mean(x), sigma = sd(x))
}))
p <- ggplot(data, aes(x = value)) +
  # shaded sigma band
  geom_rect(data = panel_stats,
            aes(xmin = mu - sigma, xmax = mu + sigma,
                ymin = 0, ymax = Inf),
            fill = "#e8a0a0", alpha = 0.4,
            inherit.aes = FALSE) +
  # histogram
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#7b9fd4", color = "white", alpha = 0.85) +
  # mu dashed line
  geom_vline(data = panel_stats,
             aes(xintercept = mu),
             linetype = "dashed", color = "#5050c0", linewidth = 0.6) +
  # labels
  geom_text(data = panel_stats,
            aes(x = mu - sigma, y = Inf, label = "-\u03c3"),
            vjust = 1.8, hjust = 0.5, color = "red", size = 5) +
  geom_text(data = panel_stats,
            aes(x = mu, y = Inf, label = "\u03bc"),
            vjust = 1.8, hjust = 0.5, color = "#5050c0", size = 5) +
  geom_text(data = panel_stats,
            aes(x = mu + sigma, y = Inf, label = "+\u03c3"),
            vjust = 1.8, hjust = 0.5, color = "red", size = 5) +
  scale_x_continuous(limits = c(-7, 9), breaks = seq(-5, 5, 5)) +
  facet_wrap(~ panel, ncol = 2) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 11),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(x = expression(EE[i]), y = "relative frequency")
p


getwd()

ggsave(
  "report/Bilder/example_grid_mu_sigma_vairation.png",
  p,
  width = 9,
  height = 6,
  dpi = 300
)