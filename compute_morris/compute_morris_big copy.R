library(sensitivity)
library(tidyverse)


# ── 0. LOAD DATA ──────────────────────────────────────────────────────────────

# load the Morris sensitivity analysis object produced by the parameter sampling
# script. This object holds the sampled parameter tuples and will later receive
# the model outputs via tell() to compute the sensitivity indices.
sa <- readRDS("compute_morris/sa_object_big.rds")

sa$X

# load the lumped statistics produced by the post-processing script.
# each row represents one annotated CWP location under one Morris parameter
# tuple, with the total detected area as the scalar model output.

lumped_statistics <- read.csv("compute_morris/lumped_stats_emme_big_2.csv")

# ── 1. NORMALISE MODEL OUTPUT ─────────────────────────────────────────────────
# Normalise the total detected area per annotated CWP location by dividing by
# the mean area across all parameter tuples for that location. This makes the
# model output comparable across locations that differ in absolute patch size —
# a location with a large CWP would otherwise dominate the sensitivity indices.

normalized_lumped_statistics <- lumped_statistics %>%
  group_by(identifier) %>%
  mutate(
    normalized_total_area = total_area / max(total_area)
  ) %>%
  ungroup()


# ── 2. COMPUTE MORRIS SENSITIVITY INDICES PER ANNOTATED CWP ──────────────────
# For each annotated CWP location, pass the normalised model output to the
# Morris object via tell() and extract the elementary effects. From the
# elementary effects, compute the three Morris sensitivity indices:
#   - mu:      mean elementary effect (sign carries direction of influence)
#   - mu.star: mean absolute elementary effect (overall parameter importance)
#   - sigma:   standard deviation of elementary effects (nonlinearity / interactions)

identifiers  <- normalized_lumped_statistics$identifier %>% unique()

lumped_statistics 




results_list <- vector("list", length(identifiers))

for (i in seq_along(identifiers)) {
  id <- identifiers[i]

  # extract the normalised model output for this CWP location, ordered by
  # Morris row index so the output aligns with the parameter tuple order
  # in the Morris object
  y <- normalized_lumped_statistics %>%
    filter(identifier == id) %>%
    arrange(morris_row_index) %>%
    pull(normalized_total_area)

  y

  # pass the model output vector to the Morris object.
  # tell() computes the elementary effects from y and the sampled parameter tuples.
  tell(sa, y)
  

  # extract the three Morris sensitivity indices from the elementary effects matrix
  mu      <- apply(sa$ee, 2, mean)
  mu.star <- apply(sa$ee, 2, function(x) mean(abs(x)))
  sigma   <- apply(sa$ee, 2, sd)


  # retrieve the class (Tributary / Non-Tributary) and mean detected area
  # for this CWP location — used for stratified plotting later
  Class     <- normalized_lumped_statistics %>% filter(identifier == id) %>% pull(Class)     %>% unique()
  mean_area <- normalized_lumped_statistics %>% filter(identifier == id) %>% pull(total_area) %>% mean(na.rm = TRUE)
  mean_deltaT <- normalized_lumped_statistics %>% filter(identifier == id) %>% pull(mean_deltaT) %>% mean(na.rm = TRUE)

  # store the sensitivity indices together with metadata for this CWP location
  results_list[[i]] <- data.frame(
    identifier = id,
    parameter  = colnames(sa$ee), # one row per parameter
    mu         = mu,
    mu.star    = mu.star,
    sigma      = sigma,
    mean_area  = mean_area,
    mean_deltaT = mean_deltaT,
    Class = Class
  )
}


# combine results from all annotated CWP locations into one table
results <- bind_rows(results_list) %>%
    mutate(
      Class = ifelse(Class == "NT", "non tributary-caused", "tributar-caused")
    )


results %>% colnames()

# ── 3. PLOT SENSITIVITY INDICES ───────────────────────────────────────────────

if (TRUE) {

# Visualise the normalised sensitivity indices stratified by CWP class
# (Tributary / Non-Tributary). Each point is one annotated CWP location.
# Locations classified as "Unshure" are excluded from the plots.

# sensitivity to slab halfwidth — stratified by CWP class
step_length_results <- results %>%
  filter(Class != "Unshure", parameter == "slab_halfwidth_m") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma)

# compute n per class based on distinct locations (for the info box)
n_per_class <- step_length_results %>%
  distinct(Class, identifier) %>%
  count(Class)

# build the annotation text for the info box
annotation_text <- paste0(
  "Morris trajectories: r = 64\n",
  paste0("N ", n_per_class$Class, " : ", n_per_class$n, collapse = "\n")
)

sensitivity_to_step_length <- ggplot(step_length_results, aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  ggtitle("Sensitivity of CWP area to Step Length") +
  annotate("label",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = annotation_text,
    size = 3
  ) +
  xlab("Morris sensitivity indices") +
  ylab("Index value") +
  scale_x_discrete(labels = c("mu.star" = expression(mu*"*"), "sigma" = expression(sigma)))

sensitivity_to_step_length

ggsave("./report/Bilder/step_length_results.png", sensitivity_to_step_length, 
       width = 6, height = 4, units = "in", dpi = 300)


# sensitivity to buffer pixel count — stratified by CWP class
buffer_px_results <- results %>%
  filter(Class != "Unshure", parameter == "buffer_px") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma)

n_per_class_buffer <- buffer_px_results %>%
  distinct(Class, identifier) %>%
  count(Class)

annotation_text_buffer <- paste0(
  "Morris trajectories: r = 64\n",
  paste0("N ", n_per_class_buffer$Class, " : ", n_per_class_buffer$n, collapse = "\n")
)

sensitivity_to_buffer_px <- ggplot(buffer_px_results, aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  ggtitle("Sensitivity of CWP area to Buffer Pixel Count") +
  annotate("label",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = annotation_text_buffer,
    size = 3
  ) +
  xlab("Morris sensitivity indices") +
  ylab("Index value") +
  scale_x_discrete(labels = c("mu.star" = expression(mu*"*"), "sigma" = expression(sigma)))

sensitivity_to_buffer_px

ggsave("./report/Bilder/buffer_px_results.png", sensitivity_to_buffer_px,
       width = 6, height = 4, units = "in", dpi = 300)


# sensitivity to temperature delta threshold — stratified by CWP class
delta_T_results <- results %>%
  filter(Class != "Unshure", parameter == "delta_T") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma)

n_per_class_delta <- delta_T_results %>%
  distinct(Class, identifier) %>%
  count(Class)

annotation_text_delta <- paste0(
  "Morris trajectories: r = 64\n",
  paste0("N ", n_per_class_delta$Class, " : ", n_per_class_delta$n, collapse = "\n")
)

sensitivity_to_delta_T <- ggplot(delta_T_results, aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  ggtitle("Sensitivity of CWP area to Temperature Delta Threshold") +
  annotate("label",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = annotation_text_delta,
    size = 3
  ) +
  xlab("Morris sensitivity indices") +
  ylab("Index value") +
  scale_x_discrete(labels = c("mu.star" = expression(mu*"*"), "sigma" = expression(sigma)))

sensitivity_to_delta_T

ggsave("./report/Bilder/delta_T_results.png", sensitivity_to_delta_T,
       width = 6, height = 4, units = "in", dpi = 300)

mu_star_vs_area <- results %>%
  filter(parameter == "slab_halfwidth_m") %>%
  ggplot(aes(x = mean_area, y = mu.star, color = Class)) +
  geom_point() +
  xlab(expression("Mean CWP area ("*m^2*")")) +
  ylab(expression(mu*"* (step length)")) +
  ggtitle(expression(mu*"* (step length) vs. mean CWP area")) +
  annotate("label",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = annotation_text,
    size = 3
  )

mu_star_vs_area

results %>% colnames()

mu_star_vs_area_deltaT <- results %>%
  filter(parameter == "slab_halfwidth_m") %>%
  ggplot(aes(x = mean_area, y = mu.star, color = mean_deltaT, shape = Class)) +
  geom_point(size = 3, stroke = 1.2) +
  scale_color_viridis_c(option = "C") +
  scale_shape_manual(values = c(16, 17)) +  # filled circle, filled triangle
  xlab(expression("Mean CWP area ("*m^2*")")) +
  ylab(expression(mu*"* (step length)")) +
  labs(color = expression(Delta*T~"(mean, "*degree*"C)"), shape = "Class") +
  ggtitle(expression(mu*"* (step length) by mean area, colored by "*Delta*"T")) +
  annotate("label",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = annotation_text,
    size = 3
  )

mu_star_vs_area_deltaT

ggsave("./report/Bilder/mu_star_vs_area_vs_deltaT.png", mu_star_vs_area_deltaT,
       width = 6, height = 4, units = "in", dpi = 300)


# fit linear model: mu.star ~ mean_area for slab_halfwidth_m
lm_data <- results %>%
  filter(parameter == "slab_halfwidth_m")


lm_fit_2 <- lm(mu.star ~ mean_area + mean_deltaT, data = lm_data)
lm_summary_2 <- summary(lm_fit_2)

lm_summary_2

} # end if (FALSE)




results %>% colnames()


mu_star_sig_of_all_plot <- results %>%
  pivot_longer(cols = mu:sigma, names_to = "score", values_to = "score_value") %>%
  filter(score != "mu") %>%
  ggplot(aes(x = parameter, y = score_value, fill = score)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.5, size = 1) +
  ylab("Index value") +
  xlab("Algorithm Parameters") +
  scale_x_discrete(limits = c("buffer_px", "slab_halfwidth_m", "delta_T"),
                    labels = c("Buffer Pixel Count", "Step Length", "Temperature Delta")) +
  scale_fill_discrete(name = "Morris Indices",
                       labels = c(expression(mu^"*"), expression(sigma)))


ggsave("./report/Bilder/mustar_sig_over_all.png",mu_star_sig_of_all_plot ,
       width = 6, height = 4, units = "in", dpi = 300)
