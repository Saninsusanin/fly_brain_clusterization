#include "community.h"

namespace kv {
    Community::Community(CommunityStorage community,
                         int communityDegree, int innerEdges, double modularity, int aggregatedNumber) {
        this->modularity = modularity;
        this->communityDegree = communityDegree;
        this->innerEdges = innerEdges;
        this->aggregationNumber = aggregatedNumber;
        this->community = std::move(community);
    }

    int Community::getCommunityId() const noexcept {
        return this->community.empty() ? 0 : (*this->community.begin())->communityId;
    }

    void Community::updateModularity(int edges) noexcept {
        this->modularity = Community::calculateCommunityModularity(edges, this->innerEdges, this->communityDegree);
    }

    void Community::updateInnerEdges(Graph& graph) noexcept {
        this->innerEdges = this->getInnerEdges(graph);
    }

    void Community::updateCommunityDegree() noexcept {
        this->communityDegree = this->getCommunityDegree();
    }

    double Community::getCommunityDensity() const noexcept {
        return this->community.size() == 1 ? 1. : double(this->innerEdges) / (0.5 * this->community.size() *
                                                                              (this->community.size() - 1.));
    }


    void Community::updateAggregationNumber() noexcept {
        this->aggregationNumber = 0;

        for (auto& currentNode : this->community) {
            this->aggregationNumber += std::dynamic_pointer_cast<LouvainNode>(currentNode)->aggregationNumber;
        }
    }


    void Community::insertNode(const pNode& node,
                               int innerEdgesFromNode,
                               int edges, int nodeAggregationNumber) {
        this->community.insert(node);
        this->communityDegree += node->getFullDegree();
        this->innerEdges += innerEdgesFromNode;
        this->aggregationNumber += nodeAggregationNumber;

        this->updateModularity(edges);
    }

    void Community::eraseNode(const pNode& node,
                              int innerEdgesFromNode,
                              int edges, int nodeAggregationNumber) {
        this->community.erase(node);
        this->innerEdges -= innerEdgesFromNode;
        this->communityDegree -= node->getFullDegree();
        this->aggregationNumber -= nodeAggregationNumber;

        this->updateModularity(edges);
    }

    int Community::getInnerEdges(const Graph &graph) const noexcept {
        int innerCommonDegree = 0;

        for (const auto& node : this->community) {
            innerCommonDegree += (node->getInnerEdgesFromNode(node->communityId, graph) + node->loops);
        }

        /// inner common degree is always an even number and equals to two times edges number
        return innerCommonDegree / 2;
    }

    int Community::getCommunityDegree() const noexcept {
        int _communityDegree = 0;

        for (const auto& node : this->community) {
            _communityDegree += node->getFullDegree();
        }

        return _communityDegree;
    }

    double Community::calculateCommunityModularity(int edgesNumber,
                                                   int _innerEdges,
                                                   int _communityDegree) noexcept {
        return (double(_innerEdges) / edgesNumber - pow(0.5 * _communityDegree / edgesNumber, 2));
    }
}