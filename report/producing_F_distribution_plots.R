library(tidyverse)
library(sensitivity)


sa <- readRDS("./compute_morris/sa_object_big.rds")
lumped_statistics <- read_csv("./compute_morris/lumped_stats_emme_big.csv")



normalized_lumped_statistics <- lumped_statistics %>%
    group_by(identifier) %>%
    mutate(
      normalized_total_area = total_area/ mean(total_area),
    ) %>% ungroup()


 
identifiers <- normalized_lumped_statistics$identifier %>% unique()




results_list <- vector("list", length(identifiers))

  y <- normalized_lumped_statistics %>%
    filter(identifier == 23) %>%
    arrange(morris_row_index) %>%
    pull(normalized_total_area)
  
  tell(sa, y)

plot(sa)

# maybe 6
ee_step_size <- sa$ee[,3]



ee_step_size <- rnorm(100,0,1) + 3




extreme_mu <- tibble(`Elementary Effect` = ee_step_size) %>%
  ggplot(aes(x = `Elementary Effect`, y = after_stat(count / sum(count)))) +
  geom_histogram(bins = 25, fill = "red", alpha = 0.5) +
  geom_vline(xintercept = mean(ee_step_size), color = "red", linetype = "dashed") +
  annotate("text", x = mean(ee_step_size), y = 0.3, 
           label = expression(mu), color = "red", hjust = -0.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Relative Frequency", 
       x = "Elementary Effect") +
  xlim(0,6)

ggsave("extreme_mu.png", , width = 10, height = 6)


ee_step_size <- rnorm(300, 0, 3)
ee_bar_step_size <- abs(ee_step_size)

difference_mu_mu_star <- tibble(
  `Elementary Effect` = c(ee_step_size, ee_bar_step_size),
  type = rep(c("EE", "|EE|"), each = 300)
) %>%
  ggplot(aes(x = `Elementary Effect`, y = after_stat(count / sum(count)), fill = type)) +
  geom_histogram(bins = 25, position = "identity", alpha = 0.5) +
  geom_vline(xintercept = mean(ee_step_size), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mean(ee_bar_step_size), color = "red", linetype = "dashed") +
  annotate("text", x = mean(ee_step_size), y = 0.1,
           label = expression(mu), color = "blue", hjust = -0.3) +
  annotate("text", x = mean(ee_bar_step_size), y = 0.1,
           label = expression(mu^"*"), color = "red", hjust = -0.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Relative Frequency",
       x = "value",
       fill = ""
      )


ggsave("report/Bilder/mu_mu_star.png", difference_mu_mu_star, width = 10, height = 6)


large_mu_large_sig <- rnorm(100,0,2) + 3
small_mu_small_sig <- rnorm(100,0,0.5)
small_mu_large_sig <- rnorm(100,0,2)
large_mu_small_sig <- rnorm(100,0,0.5) + 3



sigma_mu_comparison <- tibble(
  `Large mu large sigma` = large_mu_large_sig,
  `Small mu small sigma` = small_mu_small_sig,
  `Small mu large sigma` = small_mu_large_sig,
  `Large mu small sigma` = large_mu_small_sig
) %>% pivot_longer(
  names_to = "Combination",
  values_to = "value",
  cols = everything()
)


sigma_mu_comparison <- tibble(
  `Large Effect with nonlinearity/interactions` = large_mu_large_sig,
  `Small Effect with no nonlinearity/interactions` = small_mu_small_sig,
  `Large Effect obstructed by symetry` = small_mu_large_sig,
  `Large Effect with no nonlinearity/interactions` = large_mu_small_sig
) %>% pivot_longer(
  names_to = "Combination",
  values_to = "value",
  cols = everything()
)


stats <- sigma_mu_comparison %>%
  group_by(Combination) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  mutate(xmin = mean - sd, xmax = mean + sd,
         ypos = 0.08)     # adjust to suit your histogram scale

example_grid_of_f_distributions <- sigma_mu_comparison %>% 
  ggplot(aes(x = value, y = after_stat(count / sum(count)))) +
  geom_histogram(fill = "blue", col = "black", bins = 30, alpha = 0.5) + 
  geom_rect(data = stats,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            inherit.aes = FALSE,
            fill = "red", alpha = 0.15) +
  geom_vline(data = stats, aes(xintercept = mean), color = "blue",
             linetype = "dashed") +
  geom_text(data = stats, aes(x = mean, y = ypos, label = "μ"),
            color = "blue", hjust = -0.1) +
  geom_text(data = stats, aes(x = xmin, y = ypos, label = "-σ"),
            color = "red", hjust = 1.1) +
  geom_text(data = stats, aes(x = xmax, y = ypos, label = "+σ"),
            color = "red", hjust = -0.1) +
  facet_wrap(~ Combination , ncol = 2, axes = "all_x") +
  ylab("relative frequency")


ggsave("report/Bilder/example_grid_mu_sigma_vairation.png", example_grid_of_f_distributions, width = 11, height = 7)
