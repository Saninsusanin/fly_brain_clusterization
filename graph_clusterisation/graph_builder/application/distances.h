#ifndef APPLICATION_DISTANCES_H
#define APPLICATION_DISTANCES_H

#include <vector>
#include <iostream>
#include <cmath>
#include <functional>

#include "kernels.h"

double get_minkovskiy_distance(std::vector<double>& first_point,
                               std::vector<double>& second_point,
                               double p);
double get_minkovskiy_distance(std::vector<double>& first_point,
                               std::vector<double>& second_point,
                               std::string& max);
double get_kernel_distance(std::vector<double>& first_point,
                           std::vector<double>& second_point,
                           std::function<double(std::vector<double>&, std::vector<double>&)>& kernel);

#endif //APPLICATION_DISTANCES_H
