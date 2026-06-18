library(sf)


# move into the comparison directory, where the per-parameter-combination
# output folders are located
setwd("comparison")

# load the jaccard_similarity() helper function
source("../functions/jaccard.R")

# read the polygon outputs for parameter combination 1:
# "new" = output of the new implementation/algorithm variant
# "current" = output of the speed optimized implementation
new_combi_1     <- sf::read_sf("./param_combi_1/new/final_polys.shp")
current_combi_1 <- sf::read_sf("./param_combi_1/current/final_polys.shp")

# read the polygon outputs for parameter combination 2
new_combi_2     <- sf::read_sf("./param_combi_2/new/final_polys.shp")
current_combi_2 <- sf::read_sf("./param_combi_2/current/final_polys.shp")

# read the polygon outputs for parameter combination 3
new_combi_3     <- sf::read_sf("./param_combi_3/new/final_polys.shp")
current_combi_3 <- sf::read_sf("./param_combi_3/current/final_polys.shp")




v3 <- sf::read_sf("./param_combi_1/v3_12_core/final_polys.shp")


# compute the Jaccard similarity (intersection over union of the combined
# patch areas) between the new and current implementation, for each
# parameter combination, to assess output equivalence
jaccard_similarity(new_combi_1, current_combi_1)
jaccard_similarity(new_combi_2, current_combi_2)
jaccard_similarity(new_combi_3, current_combi_3)



jaccard_similarity(v3,current_combi_1)
jaccard_similarity(v3,current_combi_2)
jaccard_similarity(v3,current_combi_3)
