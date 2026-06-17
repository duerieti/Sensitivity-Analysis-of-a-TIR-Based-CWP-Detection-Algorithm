#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector idw_all_cells(NumericVector cell_x,        // a vector containing the x-coordinates of the pixels
                             NumericVector cell_y,        // a vector containing the y-coordinates of the pixels
                             NumericVector cell_vals,     // a vector holding the temperature values of the pixels
                             NumericMatrix idw_points,    // a matrix holding the coordinates of the reference points
                             NumericVector temperatures,  // a vector holding the temperatures of the reference points
                             double power = 2.0) {        // the decay power of the inverse distance weighting

  int ncells = cell_x.size();   // get the number of cells passed to the function
  int npts   = idw_points.nrow(); // get the number of reference points passed to the function

  NumericVector out(ncells);      // define a numeric vector with length equal to the number of cells
  LogicalVector valid_pt(npts);   // define a logical vector with length equal to the number of reference points

  for (int j = 0; j < npts; j++)                       // for each reference point temperature:
    valid_pt[j] = !std::isnan(temperatures[j]);         // store true if the temperature is valid, false if NaN

  for (int i = 0; i < ncells; i++) {                   // for every cell
    if (std::isnan(cell_vals[i])) { out[i] = NA_REAL; continue; } // if the cell temperature is NaN, paste NaN in the output vector and move to the next iteration

    double x = cell_x[i], y = cell_y[i];               // get the x and y coordinates of the current cell
    double wsum = 0.0, vwsum = 0.0;                     // initialise the cumulative sum of weights and the cumulative sum of weighted temperatures

    for (int j = 0; j < npts; j++) {                   // for every reference point
      if (!valid_pt[j]) continue;                       // if the reference point temperature was flagged as NaN, skip to the next reference point

      double dx = idw_points(j, 0) - x;                // compute the distance in the x-direction between the reference point and the cell
      double dy = idw_points(j, 1) - y;                // compute the distance in the y-direction between the reference point and the cell
      double d  = std::sqrt(dx*dx + dy*dy);             // use Pythagoras' theorem to get the Euclidean distance between the reference point and the cell

      if (d == 0.0) { vwsum = temperatures[j]; wsum = 1.0; break; } // if the distance is zero (the reference point lies directly on the cell), use its temperature directly and stop

      double w  = 1.0 / std::pow(d, power);            // compute the IDW weight as 1 / (distance ^ power)
      wsum  += w;                                       // add the weight to the cumulative sum of weights
      vwsum += w * temperatures[j];                     // add the weighted temperature to the cumulative sum of weighted temperatures
    }

    out[i] = (wsum == 0.0) ? NA_REAL : vwsum / wsum;  // set the output temperature for this cell to the weighted average; if no valid reference points were found, set to NaN
  }

  return out; // return the vector of interpolated reference temperatures for all cells
}
