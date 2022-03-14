#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include "graph_builder.h"

static void my_split(const std::string& str, std::vector<std::string>& cont, char delim = ' ')
{
    std::string token;
    std::stringstream string_stream(str);

    for (int i = 0; i < 12; ++i) {
        std::getline(string_stream, token, delim);
        cont[i] = token;
    }

}

static std::vector<std::vector<std::pair<int, std::vector<double>>>> get_hemispheres(std::string& path_to_file) {
    auto hemispheres = std::vector<std::vector<std::pair<int, std::vector<double>>>>(2);

    std::fstream csv_data(path_to_file, std::ios::in);

    if (csv_data.is_open()) {
        std::string line;
        getline(csv_data, line);
        std::vector<std::string> row(12);

        while (getline(csv_data, line)) {
            my_split(line, row, ',');

            if (row[9] == "l") {
                hemispheres[0].emplace_back(std::stoi(row[0]),
                                            std::vector<double>({std::stod(row[1]),
                                                                 std::stod(row[2]),
                                                                 std::stod(row[3])}));
            } else if (row[9] == "r") {
                hemispheres[1].emplace_back(std::stoi(row[0]),
                                            std::vector<double>({std::stod(row[1]),
                                                                 std::stod(row[2]),
                                                                 std::stod(row[3])}));
            }
        }

    } else {
        std::cout << path_to_file << std::endl;
        std::cout << "file not found" << std::endl;
    }

    return hemispheres;
}

static int get_weight(double distance, std::vector<int>& weights, double step) {
    if (weights.size() == 1) {
        return weights[0];
    }

    int left = 1;
    int right = weights.size();
    int index = ceil(right / 2.);

    while (left <= right) {

        if ((step * (index - 1) <= distance) && (distance <= step * index)) {
            return weights[weights.size() - index];
        } else if (distance > step * index) {
            left = index + 1;
        } else {
            right = index - 1;
        }

        index = ceil((right - left) / 2.) + left;
    }
}

std::vector<std::vector<int>> get_graph(std::string& path_to_file,
        double mean_size_of_cluster,
        std::vector<int>& weights,
        std::function<double(std::vector<double>&, std::vector<double>&)>& distanceFunction) {
    double distance;
    std::vector<std::vector<int>> graph;

    double step = mean_size_of_cluster / weights.size();
    auto hemispheres = get_hemispheres(path_to_file);

    for (auto& hemisphere : hemispheres) {

        for (int first_index = 0; first_index < hemisphere.size(); ++first_index) {

            for (int second_index = first_index + 1; second_index < hemisphere.size(); ++second_index) {
                distance = distanceFunction(hemisphere[first_index].second,
                        hemisphere[second_index].second);

                if (distance <= mean_size_of_cluster) {
                    graph.emplace_back(std::initializer_list<int>{hemisphere[first_index].first,
                            hemisphere[second_index].first,
                            get_weight(distance, weights, step)});
                }

            }

        }

    }

    return graph;
}

void print_graph(std::vector<std::vector<int>>& graph, std::string&path_to_file) {
    std::fstream file;
    file.open(path_to_file, std::ios::out);

    for (auto& vertex : graph) {
        file << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << std::endl;
    }
}