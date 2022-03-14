#ifndef APPLICATION_GRAPH_BUILDER_H
#define APPLICATION_GRAPH_BUILDER_H

#include <vector>

std::vector<std::vector<int>> get_graph(std::string& path_to_file,
                                        double mean_size_of_cluster,
                                        std::vector<int>& weights,
                                        std::function<double(std::vector<double>&, std::vector<double>&)>& distanceFunction);
void print_graph(std::vector<std::vector<int>>& graph, std::string& path_to_file);

#endif //APPLICATION_GRAPH_BUILDER_H
