#ifndef HGMC_LOUVAIN_ALGORITHM_H
#define HGMC_LOUVAIN_ALGORITHM_H

#include "../data_structures/nodes.h"
#include "../data_structures/community.h"
#include "../data_structures/partition.h"
#include "../data_structures/graph.h"

namespace kv {
    void Louvain(Graph& graph, bool notStochastic=true, bool notInitialized=true);
    void deterministicMoveNode(pNode& node, Graph& graph);
}

#endif //HGMC_LOUVAIN_ALGORITHM_H
