#ifndef HGMC_GENETIC_ALGORITHM_H
#define HGMC_GENETIC_ALGORITHM_H

#include "../data_structures/nodes.h"
#include "../data_structures/community.h"
#include "../data_structures/partition.h"
#include "../data_structures/graph.h"
#include "../data_structures/genetic_config.h"
#include <string>
#include <vector>
#include "louvain_algorithm.h"
#include "leiden_algorithm.h"
#include <functional>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>

namespace kv {
    class kNearestConfig {
    public:
        static const int MIN_NUMBER_OF_COMMUNITIES = 10;
        static const int MAX_NUMBER_OF_COMMUNITIES = 15;
    };


    enum Indices{
        FIRST_PARENT_ID=0,
        SECOND_PARENT_ID,
        PLACE_ID,
    };


    struct EdgeInfo {
        bool isDeleted;
        double edgeBetweenness;

        EdgeInfo() : isDeleted(false), edgeBetweenness(0) {};
        EdgeInfo(bool isDeleted, int edgeBetweenness) : isDeleted(isDeleted), edgeBetweenness(edgeBetweenness) {};
    };


    struct NodeInfo {
        bool isUsed;
        int numberOfPaths;
        int level;
        double ebStorage;
        int vertexId;
    };


    void kNearestNeighboursPartition(Graph& graph);
    void geneticAlgorithm(Graph& graph, std::function<void(Graph&)>& method,
                          std::vector<std::function<void(kv::Graph&, GeneticConfig&)>>& mutations,
                          GeneticConfig& config);
    void randomInitialization(Graph& graph, int numberOfClusters);
    std::vector<std::vector<int>> GirvanNewman(Graph& graph, std::vector<pNode>& community) noexcept;
    std::vector<std::vector<int>> randomSplit(std::vector<pNode>& community) noexcept;

    void mergingMutate(Graph& graph, GeneticConfig& config) noexcept;
    void extractingMutate(Graph& graph, GeneticConfig& config) noexcept;
    void separatingMutate(Graph& graph, GeneticConfig& config,
                          std::function<std::vector<std::vector<int>>(std::vector<kv::pNode>&)>& engine) noexcept;
    pGraph cross(const Graph& firstGraph, const Graph& secondGraph) noexcept;
}

#endif //HGMC_GENETIC_ALGORITHM_H
