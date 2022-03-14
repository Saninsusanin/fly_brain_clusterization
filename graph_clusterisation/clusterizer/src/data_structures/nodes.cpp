#include "nodes.h"

#include <utility>
#include <memory>
#include "graph.h"

namespace kv {
    Node::Node(int vertexId, int communityId, int degree, int loops) {
        this->vertexId = vertexId;
        this->communityId = communityId;
        this->degree = degree;
        this->loops = loops;
    }

    int Node::getFullDegree() const noexcept {
        return this->degree + 2 * this->loops;
    }

    int Node::getInnerEdgesFromNode(int _communityId, const Graph &graph) const noexcept {
        int innerNodeDegree = 0;

        for (const auto& neighbourInfo : (*graph.adjacencyList)[this->vertexId - 1]) {
            const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
            if (_communityId == neighbour->communityId) {
                innerNodeDegree += neighbourInfo.edges;
            }
        }

        return innerNodeDegree + this->loops;
    }

    std::unordered_map<int, int> Node::getInnerEdgesFromNode(const Graph &graph) const noexcept {
        std::unordered_map<int, int> innerEdges{};
        for (const auto& neighbourInfo : (*graph.adjacencyList)[this->vertexId - 1]) {
            const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
            int _communityId = neighbour->communityId;
            if (innerEdges.find(_communityId) == innerEdges.end()) {
                innerEdges[_communityId] = this->loops + neighbourInfo.edges;
            } else {
                innerEdges[_communityId] += neighbourInfo.edges;
            }
        }
        return innerEdges;
    }

    pNode Node::copy() const noexcept {
        return std::make_shared<Node>(this->vertexId, 0, this->degree, this->loops);
    }


    LouvainNode::LouvainNode(int vertexId, int communityId, int degree, int loops,
                             std::vector<pNode> aggregatedVertices, int aggregatedNumber)
            : Node(vertexId, communityId, degree, loops) {
        this->aggregationNumber = aggregatedNumber;
        this->aggregatedVertices = std::move(aggregatedVertices);
    }

    LouvainNode::LouvainNode(const pNode& node) : Node(node->vertexId, node->communityId, node->degree, node->loops) {
        this->aggregationNumber = 1;
        this->aggregatedVertices = std::move(std::vector<pNode>());
    }

    pNode LouvainNode::copy() const noexcept {
        return std::make_shared<LouvainNode>(this->vertexId, 0, this->degree, this->loops,
                                             std::vector<pNode>(), this->aggregationNumber);
    }
}
