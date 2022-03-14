#ifndef APPLICATION_KERNELS_H
#define APPLICATION_KERNELS_H

#include <vector>
#include <iostream>
#include <cmath>

double polynomial_kernel(std::vector<double>& first_point,
                         std::vector<double>& second_point,
                         double gamma, double bias, double power);
double sigmoid_kernel(std::vector<double>& first_point,
                      std::vector<double>& second_point,
                      double gamma, double bias);
double exponential_kernel(std::vector<double>& first_point,
                          std::vector<double>& second_point,
                          double gamma, double p);

#endif //APPLICATION_KERNELS_H
