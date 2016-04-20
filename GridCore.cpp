#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix GridCore(NumericMatrix image, int g, int i, int j) {
   return image(Range(i*g-g,i*g-1),Range(j*g-g,j*g-1));
}
