#include "kernels.h"
#include "distances.h"

double polynomial_kernel(std::vector<double>& first_point,
        std::vector<double>& second_point,
        double gamma, double bias, double power) {
    double value = 0;

    if (first_point.size() != second_point.size()) {
        std::cerr << "vectors have different number of coordinates" << std::endl;
    }

    for (int index = 0; index < first_point.size(); ++index) {
        value += first_point[index] * second_point[index];
    }

    return pow(gamma * value + bias, power);
}

double sigmoid_kernel(std::vector<double>& first_point,
                      std::vector<double>& second_point,
                      double gamma, double bias) {
    double value = 0;

    if (first_point.size() != second_point.size()) {
        std::cerr << "vectors have different number of coordinates" << std::endl;
    }

    for (int index = 0; index < first_point.size(); ++index) {
        value += first_point[index] * second_point[index];
    }

    return tanh(gamma * value + bias);
}

double exponential_kernel(std::vector<double>& first_point,
        std::vector<double>& second_point,
        double gamma, double p) {
    return exp(-gamma * get_minkovskiy_distance(first_point, second_point, p));
}