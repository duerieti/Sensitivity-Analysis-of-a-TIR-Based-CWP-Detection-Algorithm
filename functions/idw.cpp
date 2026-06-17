#include <Rcpp.h>
#include <cmath>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector idw_all_cells(NumericVector cell_x, // a vector containing the x-coordinates of the pixels
                            NumericVector cell_y, // a vector containing the y-coordinates of the pixels
                            NumericVector cell_vals, // a vector holding the temperature values of the pixels
                            NumericMatrix idw_points, // a matrix holding the coordinates of the reference points
                            NumericVector temperatures, // a vector holding the temperatures of the reference points
                            double power   = 2.0, // the decay power of the inverse distance weighting
                            int    threads = 4) { // number of threads for parallelism
  
  int ncells = cell_x.size(); // get the number of cells passed to the function
  int npts   = idw_points.nrow(); // get the number of reference points passed to the function
  NumericVector out(ncells); // define a nummeric vector with length of number of cells

  LogicalVector valid_pt(npts); // define a vector with length of equal to the number of points
  for (int j = 0; j < npts; j++) // for each temperature of the reference points:
    valid_pt[j] = !std::isnan(temperatures[j]); // if the temperature is NaN store a true in the vector else false

  #pragma omp parallel for num_threads(threads) schedule(static) // parallelism for this loop
  for (int i = 0; i < ncells; i++) {  // for every cell 
    if (std::isnan(cell_vals[i])) { out[i] = NA_REAL; continue; } // if the temperature value of the cell is NaN, paste NaN in the output vector and on to the next itteration
    double x = cell_x[i], y = cell_y[i]; // else get the x and the y coordinates of the cells from the vectors
    double wsum = 0.0, vwsum = 0.0; // initiate two variables where the sum of the weights and the weighted sum of the temperatures will be collected
    for (int j = 0; j < npts; j++) { // for every reference point
      if (!valid_pt[j]) continue; // if the temperature was flaged as NaN, then continue to the next itteration
      double dx = idw_points(j, 0) - x; // else compute the distance in x-direction between reference point and cell
      double dy = idw_points(j, 1) - y; // and the distance in the y-direction between reference point and cell
      double d  = std::sqrt(dx*dx + dy*dy); // then use pythagoras theorem to get the distance between reference point and cell
      if (d == 0.0) { vwsum = temperatures[j]; wsum = 1.0; break; } // if the distance is zero (meaning, point is dirreclty on the cell)
      double w  = 1.0 / std::pow(d, power); // compute the weight as 1/(distance)^power
      wsum  += w; vwsum += w * temperatures[j]; // ad the weight to the cumulative sum of weights, add the weighted temperature to the cumulative sum of weighted temperatures
    }
    out[i] = (wsum == 0.0) ? NA_REAL : vwsum / wsum; // set the final output temperature for this pixel position to vwsum / wsum (which is the weighted temperature)
  }
  return out; // return the computed weighted temperatures for the cells that were passed to the function
}
