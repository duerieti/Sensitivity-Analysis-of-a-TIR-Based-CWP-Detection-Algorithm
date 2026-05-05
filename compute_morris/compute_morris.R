library(sensitivity)
library(tidyverse)


getwd()

sa <- readRDS("morris_screening/sa_object.rds")
lumped_statistics <- read.csv("compute_morris/lumped_stats_emme.csv")

lumped_statistics %>% colnames()


normalized_lumped_statistics <- lumped_statistics %>%
    group_by(identifier) %>%
    mutate(
      normalized_total_area = total_area/ mean(total_area),
      normalized_temperature = total_averaged_temperature / mean(total_averaged_temperature )
    ) %>% ungroup()



identifiers <- normalized_lumped_statistics$identifier %>% unique()

results_list <- vector("list", length(identifiers))

for (i in seq_along(identifiers)) {
  id <- identifiers[i]
  
  y <- normalized_lumped_statistics %>%
    filter(identifier == id) %>%
    arrange(morris_row_index) %>%
    pull(normalized_temperature)
  
  tell(sa, y)
  
  mu      <- apply(sa$ee, 2, mean)
  mu.star <- apply(sa$ee, 2, function(x) mean(abs(x)))
  sigma   <- apply(sa$ee, 2, sd)
  
  Class      = normalized_lumped_statistics %>%  filter(identifier == id) %>% pull(Class) %>% unique()
  mean_area  = normalized_lumped_statistics %>%  filter(identifier == id) %>% pull(total_area) %>% median(na.rm = TRUE)
  mean_temp  = normalized_lumped_statistics %>%  filter(identifier == id) %>% pull(total_averaged_temperature) %>% median(na.rm = TRUE)

  results_list[[i]] <- data.frame(
    identifier = id,
    parameter  = colnames(sa$ee),
    mu         = mu,
    mu.star    = mu.star,
    sigma      = sigma,
    Class = Class,
    mean_area = mean_area,
    mean_temp = mean_temp
   
  )
}


results <- bind_rows(results_list)

results_normalized <- results %>% as.tibble() %>% drop_na() %>%
  group_by(identifier) %>%
  mutate(
    normalized_mu = mu / sum(mu),
    normalized_mu.star = mu.star / sum(mu.star),
    normalized_sigma = sigma / sum(sigma)

)


results_normalized  %>% as.tibble() %>%
  filter(Class != "Unshure", parameter == "slab_halfwidth_m") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = normalized_mu:normalized_sigma) %>%
  ggplot(aes(y = value, x = quantity, fill = Class)) + geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))


results_normalized  %>% as.tibble() %>%
  filter(Class != "Unshure", parameter == "buffer_px") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = normalized_mu:normalized_sigma) %>%
  ggplot(aes(y = value, x = quantity, fill = Class)) + geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))




results_normalized %>%
  filter(parameter == "buffer_px") %>%
  ggplot(aes(y = mu.star, x = mean_area, color = Class)) + geom_point()
