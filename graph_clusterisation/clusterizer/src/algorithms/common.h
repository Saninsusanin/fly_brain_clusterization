#ifndef HGMC_COMMON_H
#define HGMC_COMMON_H

#include "../data_structures/nodes.h"
#include "../data_structures/community.h"
#include "../data_structures/graph.h"
#include "../data_structures/partition.h"

namespace kv {
    double getInitialDelta(std::vector<pCommunity>& partition, const pNode& node,
                           bool isAlreadyInSingleton, int edges,
                           std::unordered_map<int, int>& innerEdges);
    double getCurrentDelta(const pCommunity& community, const pNode& node,
                           int edges, std::unordered_map<int, int>& innerEdges);
    double getCurrentSingletonDelta(const pNode& node, int edges);
    void swapCommunities(Graph& graph, pNode& currentSingletonCommunity);

    void aggregateNodes(Graph& graph);
    std::vector<pCommunity> getSingletonPartition(Graph& graph);
    std::vector<pCommunity> flat(const Partition& partition);

}

#endif //HGMC_COMMON_H
