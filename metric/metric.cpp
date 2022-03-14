#include "metric.h"


kv::partition read_partition(std::string& path_to_partition) {
    int community_id = 1;
    std::string line;
    kv::partition partition;
    std::fstream partition_stream(path_to_partition, std::ios::in);

    if (partition_stream.is_open()) {
        ///get rid of first line with metadata
        std::getline(partition_stream, line, '\n');

        while (std::getline(partition_stream, line, '\n')) {
            std::string token;
            std::stringstream string_stream(line);
            partition[community_id] = std::unordered_map<int, int>();

            while (std::getline(string_stream, token, ' ')) {
                partition[community_id][std::stoi(token)] = std::stoi(token);
            }

            ++community_id;
        }
    } else {
        std::cout << "file not found" << std::endl;
    }

    return partition;
}


double kv_metric(kv::partition& real_partition, kv::partition& current_partition) {
    ///get amount of points in brain
    int number_of_vertices = 0;

    for (auto& community : real_partition) {
        number_of_vertices += community.second.size();
    }

    double metric_value = 0;

    for (auto& community_from_current : current_partition) {
        /**
            * key - cluster id from real_partition
            * value - number of vertices in this cluster
        * **/
        std::unordered_map<int, int> current_values;

        ///initialize current_values
        for (auto& community_from_real : real_partition) {
            current_values.insert({community_from_real.first, 0});
        }

        for (auto& vertex : community_from_current.second) {
            for (auto& community_from_real : real_partition) {

                if (community_from_real.second.find(vertex.first) != community_from_real.second.end()) {
                    current_values[community_from_real.first] += 1;
                }

            }
        }

        ///getting of maximum size of intersection
        int maximum_amount_of_vertices = 0;
        int res = 0;

        for (auto& item : current_values) {
            res += item.second;
            maximum_amount_of_vertices = maximum_amount_of_vertices < item.second ?
                    item.second :
                    maximum_amount_of_vertices;
        }

        metric_value += maximum_amount_of_vertices;
    }

    return metric_value / number_of_vertices;
}