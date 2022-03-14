#ifndef METRIC_METRIC_H
#define METRIC_METRIC_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>

namespace kv{
    using partition = std::unordered_map<int, std::unordered_map<int, int>>;
}

kv::partition read_partition(std::string& path_to_partition);
double kv_metric(kv::partition& real_partition, kv::partition& current_partition);

#endif //METRIC_METRIC_H
