library(sensitivity)

sa <- morris(
  model   = NULL,                    # NULL = just generate the sample
  factors = c("buffer_width", "step_length"),
  r       = 100,                     # trajectories
  design  = list(
    type      = "oat",
    levels    = 10,
    grid.jump = 5
  ),
  binf = c(100, 10),                 # lower bounds
  bsup = c(1000, 100)                # upper bounds
)

# Your parameter combinations — ready for GNU Parallel
write.table(
  cbind(index = 1:nrow(sa$X), sa$X),
  "params.txt",
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(sa, "sa_object.rds")
