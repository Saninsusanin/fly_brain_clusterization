#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>

enum coords {
    X, Y, Z
};

using point_t = std::array<double, 3>;
using point_process_t = std::vector<point_t>;


static double get_euclidean_distance(point_t& first_point,
                                     point_t& second_point) {
    return sqrt((first_point[coords::X] - second_point[coords::X]) * \
    (first_point[coords::X] - second_point[coords::X]) + \
    (first_point[coords::Y] - second_point[coords::Y]) * \
    (first_point[coords::Y] - second_point[coords::Y]) + \
    (first_point[coords::Z] - second_point[coords::Z]) * \
    (first_point[coords::Z] - second_point[coords::Z]));
}


static std::shared_ptr<point_t> get_knn(point_t& x, point_process_t& process, int k) {
    std::vector<std::pair<double, std::shared_ptr<point_t>>> distance_storage;
    distance_storage.reserve(process.size());
    distance_storage.resize(distance_storage.capacity());

    for (int i = 0; i < process.size(); ++i) {
        distance_storage[i] = {get_euclidean_distance(x, process[i]), std::make_shared<point_t>(process[i])};
    }

    std::sort(distance_storage.begin(), distance_storage.end(),
            [](std::pair<double, std::shared_ptr<point_t>>& a,
                    std::pair<double, std::shared_ptr<point_t>>& b) {
        return a.first < b.first;
    });

    return distance_storage[k].second;
}


double lambda_estimate(point_t& u, point_process_t& process, int k) {
    double lambda = 0;

    for (auto& point : process) {
        lambda += get_euclidean_distance(u, *get_knn(point, process, k)) ;
    }

    return (int(process.size()) * k - 1) / ((4 * M_PI / 3) * lambda);
}


bool point_criterion(point_process_t& process_1, int k_1,
        point_process_t& process_2, int k_2,
        double significance_level,
        point_t& u) {
    double lambda_1 = lambda_estimate(u, process_1, k_1);
    double lambda_2 = lambda_estimate(u, process_2, k_2);
    double R_12 = lambda_2 / lambda_1;

    return (std::pow(R_12, significance_level / 2) <= R_12) and \
    (R_12 <= std::pow(R_12, 1 - significance_level / 2));
}


bool process_criterion(point_process_t& process_1, int k_1,
        point_process_t& process_2, int k_2,
        double significance_level, double threshold) {
    double score = 0;

    for (auto& point : process_1) {
        score += point_criterion(process_1, k_1, process_2, k_2, significance_level, point);
    }

    return threshold <= (score / int(process_1.size()));
}


void my_split(std::string& str, std::vector<std::string>& cont, char delim = ' ')
{
    size_t pos;
    std::string token;
    std::stringstream string_stream(str);

    while ((pos = str.find(delim)) != std::string::npos) {
        token = str.substr(0, pos);
        cont.push_back(token);
        str.erase(0, pos + 1);
    }

}


std::unordered_map<int, point_t> get_point_data(std::string& path_to_file) {
    auto point_data = std::unordered_map<int, point_t>();

    std::fstream csv_data(path_to_file, std::ios::in);

    if (csv_data.is_open()) {
        std::string line;
        getline(csv_data, line);

        while (getline(csv_data, line)) {
            std::vector<std::string> row;
            my_split(line, row, ',');
            point_data.insert({std::stoi(row[0]),
                         {std::stod(row[1]),
                          std::stod(row[2]),
                          std::stod(row[3])}});
        }

    } else {
        std::cout << path_to_file << std::endl;
        std::cout << "file not found" << std::endl;
    }

    return point_data;
}


std::unordered_map<int, std::vector<int>> get_cluster_data(std::string& path_to_file) {
    int cluster_id = 0;
    auto cluster_data = std::unordered_map<int, std::vector<int>>();

    std::fstream partition_data(path_to_file, std::ios::in);

    if (partition_data.is_open()) {
        std::string line;
        getline(partition_data, line);

        while (getline(partition_data, line)) {
            std::vector<std::string> row;
            my_split(line, row);
            std::vector<int> vertices;
            vertices.reserve(row.size());

            for (auto& element : row) {
                vertices.push_back(std::stoi(element));
            }

            cluster_data.insert({cluster_id++, vertices});
        }

    } else {
        std::cout << path_to_file << std::endl;
        std::cout << "file not found" << std::endl;
    }

    return cluster_data;
}


std::vector<point_process_t> clusters_to_processes(std::vector<std::unordered_map<int, point_t>>& points,
        std::vector<std::unordered_map<int, std::vector<int>>>& clusters) {
    std::vector<point_process_t> processes;
    processes.reserve(clusters[0].size());
    processes.resize(processes.capacity());

    for (int data_id = 0; data_id < points.size(); ++data_id) {
        for (int cluster_id = 0; cluster_id < processes.size(); ++cluster_id) {
            processes[cluster_id].reserve(clusters[data_id][cluster_id].size());

            for (auto point_id : clusters[data_id][cluster_id]) {
                processes[cluster_id].push_back(points[data_id][point_id]);
            }
        }
    }

    return processes;
}


std::vector<std::vector<int>> test(std::vector<point_process_t>& first, std::vector<point_process_t>& second) {
    std::vector<std::vector<int>> matrix;
    matrix.reserve(first.size());
    matrix.resize(matrix.capacity());


    for (int i = 0; i < first.size(); ++i) {
        matrix[i].reserve(second.size());

        for (int j = 0; j < second.size(); ++j) {
            matrix[i].push_back(int(process_criterion(first[i], 9, second[j], 9, 0.05, 0.5)));
        }
    }

    return matrix;
}


int main() {
    auto path_to_csv = std::string("../pipeline/zonetab-mel-animal-2-l/data.csv");
    auto point_data = get_point_data(path_to_csv);
    auto path_to_genetic_partition = std::string("../pipeline/zonetab-mel-animal-2-l/genetic_partition.out");
    auto genetic_partition = get_cluster_data(path_to_genetic_partition);
    auto path_to_real_partition = std::string("../pipeline/zonetab-mel-animal-2-l/real_partition.out");
    auto real_partition = get_cluster_data(path_to_real_partition);
    std::vector<std::unordered_map<int, point_t>> points = {point_data};
    std::vector<std::unordered_map<int, std::vector<int>>> clusters = {real_partition};
    auto mel_2 = clusters_to_processes(points, clusters);

    path_to_csv = std::string("../pipeline/zonetab-mel-animal-3-l/data.csv");
    point_data = get_point_data(path_to_csv);
    path_to_genetic_partition = std::string("../pipeline/zonetab-mel-animal-3-l/genetic_partition.out");
    genetic_partition = get_cluster_data(path_to_genetic_partition);
    path_to_real_partition = std::string("../pipeline/zonetab-mel-animal-3-l/real_partition.out");
    real_partition = get_cluster_data(path_to_real_partition);
    points = {point_data};
    clusters = {real_partition};
    auto mel_3 = clusters_to_processes(points, clusters);

    auto matrix = test(mel_2, mel_3);

    return 0;
}