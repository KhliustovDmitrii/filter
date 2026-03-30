#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

namespace filter::math
{
void invert_matrix(std::vector<double> &dest, std::vector<double> &source, int dim);
void matrix_square_root(std::vector<double> &dest, std::vector<double> &source, int dim);
}; // filter::math
#endif