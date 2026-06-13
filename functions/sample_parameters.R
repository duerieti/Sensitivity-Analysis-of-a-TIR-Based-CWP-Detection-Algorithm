library(sensitivity) # support for a range of sensitivity analysis methods

# sample parameter tupples from a morris design 
morrisDesign <- morris(
    model   = NULL, # passing null to the model slot tells the library that there is no model to set up.
    factors = c("buffer_px", "step_length"), # the two parameters subjectto variation
    r       = 44, # 44 trajectories in sampling (the optimized sampling strategy for morris deploys trajectories)
    design  = list(type = "oat", levels = 6, grid.jump = 3), 
    # the grid has 6 levels (in normalized space {0, 0.2, 0.4, 0.6, 0.8, 1.0} )
    # to compute the elementary effects EE_i (partial derivatives), 3 grid points will be jumped.
    # this corresponds to a Delta 3/5 in the normalized space.
  
    binf    = c(2, 400), # the lower bounds of the grid in x and y direction (2D grid)
    bsup    = c(7, 1200) # the upper bounds of the grid in x and y direction (2D grid)
)


# write the parameter tupples to a table and write it as params.txt
write.table(
  cbind(index = 1:nrow(morrisDesign$X), morrisDesign$X),
  "params.txt",
  row.names = FALSE,
  col.names = FALSE
)

# also save the sensitivity analysis object as it is later
# needed to compute the morris sensitivity indices.

saveRDS(morrisDesign, "sa_object.rds")
