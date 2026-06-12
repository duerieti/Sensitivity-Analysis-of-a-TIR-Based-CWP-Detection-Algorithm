library(tidyverse)

current_algo_execution_speed <- c(7394.236, 5046.987, 4620.728)
new_algo_execution_speed     <- c(2010.463, 1992.954, 2142.211)
step_lengths                 <- c(200, 800, 1200)

# combine into a long-format tibble, converting seconds to minutes
exec_times <- tibble(
  step_length = rep(step_lengths, 2),
  runtime_min = c(current_algo_execution_speed, new_algo_execution_speed) / 3600,
  algorithm   = rep(c("v1 (current)", "v2 (speed optimized)"), each = length(step_lengths))
)

p <- ggplot(exec_times, aes(x = step_length, y = runtime_min, color = algorithm)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    x = "Segment length (meters)",
    y = "Execution time (hours)",
    color = "Implementation"
  ) +
  theme_minimal() + ylim(0,2.2)+ xlim(200, 1250)

p + scale_x_continuous(breaks = seq(200, 1200, by = 200))
