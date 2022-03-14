#ifndef HGMC_NODES_H
#define HGMC_NODES_H

#include <set>
#include <queue>
#include <vector>
#include <random>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <fstream>

namespace kv {
    struct Node;
    struct LouvainNode;
    struct Graph;

    using pNode = std::shared_ptr<Node>;
    using pLouvainNode = std::shared_ptr<LouvainNode>;

    struct Node {
        int vertexId;
        int communityId;
        int degree;
        int loops;

        explicit Node(int vertexId=0, int communityId=0, int degree=0, int loops=0);
        virtual ~Node() = default;

        int getFullDegree() const noexcept;
        int getInnerEdgesFromNode(int _communityId, const Graph& graph) const noexcept;
        std::unordered_map<int, int> getInnerEdgesFromNode(const Graph& graph) const noexcept;

        virtual pNode copy() const noexcept;
    };

    struct LouvainNode : public Node {
        int aggregationNumber;
        std::vector<pNode> aggregatedVertices;

        explicit LouvainNode(int vertexId=0, int communityId=0, int degree=0, int loops=0,
                             std::vector<pNode> aggregatedVertices=std::vector<pNode>(),
                             int aggregationNumber=0);
        explicit LouvainNode(const pNode& node);

        pNode copy() const noexcept final;
    };
}

#endif //HGMC_NODES_H