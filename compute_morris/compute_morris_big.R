library(sensitivity)
library(tidyverse)


# ── 0. LOAD DATA ──────────────────────────────────────────────────────────────

# load the Morris sensitivity analysis object produced by the parameter sampling
# script. This object holds the sampled parameter tuples and will later receive
# the model outputs via tell() to compute the sensitivity indices.
sa <- readRDS("compute_morris/sa_object_big.rds")

# load the lumped statistics produced by the post-processing script.
# each row represents one annotated CWP location under one Morris parameter
# tuple, with the total detected area as the scalar model output.
lumped_statistics <- read.csv("compute_morris/lumped_stats_emme_big.csv")


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

  # store the sensitivity indices together with metadata for this CWP location
  results_list[[i]] <- data.frame(
    identifier = id,
    parameter  = colnames(sa$ee), # one row per parameter
    mu         = mu,
    mu.star    = mu.star,
    sigma      = sigma,
    Class      = Class,
    mean_area  = mean_area
  )
}

# combine results from all annotated CWP locations into one table
results <- bind_rows(results_list)


# ── 3. NORMALISE SENSITIVITY INDICES ─────────────────────────────────────────
# Normalise mu.star and sigma within each CWP location by dividing by their
# sum across parameters. This converts absolute sensitivity indices into
# relative parameter contributions (shares summing to 1), making it possible
# to compare the relative importance of parameters across CWP locations that
# may have very different absolute sensitivities.



# ── 4. PLOT SENSITIVITY INDICES ───────────────────────────────────────────────
# Visualise the normalised sensitivity indices stratified by CWP class
# (Tributary / Non-Tributary). Each point is one annotated CWP location.
# Locations classified as "Unshure" are excluded from the plots.

# sensitivity to slab halfwidth — stratified by CWP class
step_length_results <- results %>%
  filter(Class != "Unshure", parameter == "slab_halfwidth_m") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma) %>%
  ggplot(aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  ggtitle("Sensitivity to Slab Halfwidth")


step_length_results

ggsave("step_length_results.png", step_length_results)

# sensitivity to buffer pixel count — stratified by CWP class
results %>%
  filter(Class != "Unshure", parameter == "buffer_px") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma) %>%
  ggplot(aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

# sensitivity to temperature delta threshold — stratified by CWP class
results %>%
  filter(Class != "Unshure", parameter == "delta_T") %>%
  pivot_longer(names_to = "quantity", values_to = "value", cols = mu.star:sigma) %>%
  ggplot(aes(y = value, x = quantity, fill = Class)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))

# mu.star across all annotated CWP locations and parameters — gives an overall
# picture of which parameters drive the most variation in detected CWP area
mu_star_res <- results_normalized %>%
  pivot_longer(names_to = "measure", values_to = "value", cols = mu:sigma) %>%
  filter(measure == "mu.star") %>%
  ggplot(aes(x = parameter, y = value)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  ggtitle("mu.star across all annotated CWP locations")


mu_star_res
ggsave("mu_star_res.png", mu_star_res)
