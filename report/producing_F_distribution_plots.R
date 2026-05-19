library(tidyverse)
library(sensitivity)

getwd()

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



ee_step_size <- rnorm(100,0,1) - 3



extreme_mu <- tibble(`Elementary Effect` = ee_step_size) %>%
  ggplot(aes(x = `Elementary Effect`, y = after_stat(count / sum(count)))) +
  geom_histogram(bins = 25) +
  geom_vline(xintercept = mean(ee_step_size), color = "red", linetype = "dashed") +
  annotate("text", x = mean(ee_step_size), y = 0.3, 
           label = expression(mu), color = "red", hjust = -0.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Relative Frequency", 
       x = "Elementary Effect",
       title = expression("Distribution F an Elementary Effect approximated by Histogram ")) +
       xlim(-6, 0)


ee_step_size <- rnorm(300, 0, 3)
ee_bar_step_size <- abs(ee_step_size)

tibble(
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
       x = "Elementary Effect",
       fill = "",
       title = "Distribution F and G of an elementary Effect")


