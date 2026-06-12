# Compute the Jaccard similarity (intersection over union) between two sets
# of polygons, based on their combined (dissolved) spatial extent.
#
# Args:
#   poly_a, poly_b: sf objects (or sfc geometries) containing polygons to compare
#   make_valid:     if TRUE (default), repair invalid geometries before
#                    computing union/intersection
#
# Returns:
#   A single numeric value in [0, 1]. Returns 0 if the union area is zero
#   (i.e., both inputs are empty).

jaccard_similarity <- function(poly_a, poly_b, make_valid = TRUE) {

  # dissolve each set into a single (multi)polygon
  union_a <- sf::st_union(poly_a)
  union_b <- sf::st_union(poly_b)

  if (make_valid) {
    union_a <- sf::st_make_valid(union_a)
    union_b <- sf::st_make_valid(union_b)
  }

  # compute intersection and union geometries
  intersection_geom <- sf::st_intersection(union_a, union_b)
  union_geom        <- sf::st_union(union_a, union_b)

  # compute areas
  intersection_area <- sum(as.numeric(sf::st_area(intersection_geom)))
  union_area        <- sum(as.numeric(sf::st_area(union_geom)))

  # avoid division by zero if both inputs are empty
  if (union_area == 0) return(0)

  intersection_area / union_area
}

