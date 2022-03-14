#include "distances.h"

double get_minkovskiy_distance(std::vector<double>& first_point,
                              std::vector<double>& second_point,
                              double p) {
    double dist = 0;

    if (first_point.size() != second_point.size()) {
        std::cerr << "vectors have different number of coordinates" << std::endl;
    }

    for (int index = 0; index < first_point.size(); ++index) {
        dist += pow(fabs(first_point[index] - second_point[index]), p);
    }

    return pow(dist, 1/p);
}

double get_minkovskiy_distance(std::vector<double>& first_point,
                               std::vector<double>& second_point,
                               std::string& max) {
    double dist = 0;

    if (first_point.size() != second_point.size()) {
        std::cerr << "vectors have different number of coordinates" << std::endl;
    }

    for (int index = 0; index < first_point.size(); ++index) {
        double tmp = fabs(first_point[index] - second_point[index]);

        if (dist < tmp) {
            dist = tmp;
        }

    }

    return dist;
}

double get_kernel_distance(std::vector<double>& first_point,
                           std::vector<double>& second_point,
                           std::function<double(std::vector<double>&, std::vector<double>&)>& kernel) {
    std::vector<double> tmp;
    tmp.reserve(first_point.size());

    if (first_point.size() != second_point.size()) {
        std::cerr << "vectors have different number of coordinates" << std::endl;
    }

    for (int index = 0; index < first_point.size(); ++index) {
        tmp[index] = first_point[index] - second_point[index];
    }

    return sqrt(kernel(tmp, tmp));
}
