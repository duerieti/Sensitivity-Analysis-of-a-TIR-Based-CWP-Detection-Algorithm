library(sensitivity)

morrisDesign <- morris(
    model   = NULL,
    factors = c("buffer_px", "slab_halfwidth_m"),
    r       = 66,
    design  = list(type = "oat", levels = 6, grid.jump = 3),
    binf    = c(2, 400),
    bsup    = c(7, 1200)
)


# Your parameter combinations — ready for GNU Parallel
write.table(
  cbind(index = 1:nrow(morrisDesign$X), morrisDesign$X),
  "params.txt",
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(morrisDesign, "sa_object.rds")
