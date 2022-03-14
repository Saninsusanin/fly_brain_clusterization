#ifndef HGMC_LEIDEN_ALGORITHM_H
#define HGMC_LEIDEN_ALGORITHM_H

#include "../data_structures/nodes.h"
#include "../data_structures/community.h"
#include "../data_structures/partition.h"
#include "../data_structures/graph.h"

namespace kv {
    class LeidenConfig {
    public:
        constexpr static double TETHA = 0.05;
    };

    void Leiden(Graph& graph);
    void moveNodesFast(Graph& graph);
    void mergeNodesSubset(Graph& graph, const Community& subset);
    void refinePartition(Graph& refineGraph, const Graph& moveGraph);
    void updateMovePartition(const Graph& refineGraph, Graph& moveGraph);
}

#endif //HGMC_LEIDEN_ALGORITHM_H
