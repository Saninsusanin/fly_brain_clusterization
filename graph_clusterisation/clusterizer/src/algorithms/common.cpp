#include "../data_structures/nodes.h"
#include "common.h"

#include <unordered_map>

namespace kv {

    double getInitialDelta(std::vector<pCommunity>& partition, const pNode& node,
                           bool isAlreadyInSingleton, int edges,
                           std::unordered_map<int, int>& innerEdges) {
        if (isAlreadyInSingleton) {
            double initialDelta = -partition.back()->modularity;
            partition.pop_back();
            return initialDelta;
        }

        auto nodeCommunity = partition[node->communityId - 1];
        double oldModularity = nodeCommunity->modularity;

        nodeCommunity->eraseNode(node, innerEdges[node->communityId],
                                 edges, std::dynamic_pointer_cast<LouvainNode>(node)->aggregationNumber);

        return nodeCommunity->modularity - oldModularity;
    }


    double getCurrentDelta(const pCommunity& community, const pNode& node,
                           int edges, std::unordered_map<int, int>& innerEdges) {
        double oldModularity = community->modularity;

        ///recalculate value of metric on current community
        double newModularity = Community::calculateCommunityModularity(edges,
                                                                       community->innerEdges + innerEdges[community->getCommunityId()],
                                                                       community->communityDegree + node->getFullDegree());

        return newModularity - oldModularity;
    }


    double getCurrentSingletonDelta(const pNode& node, int edges) {
        return Community::calculateCommunityModularity(edges, node->loops, node->getFullDegree());
    }


///swap last community with current singleton community
    void swapCommunities(Graph& graph, pNode& currentSingletonCommunity) {

        ///changing of community id of the last community in partition to id of currentSingletonCommunity
        for (auto& currentNode : graph.partition->partition.back()->community) {
            currentNode->communityId = currentSingletonCommunity->communityId;
        }

        ///swap
        std::swap(graph.partition->partition[currentSingletonCommunity->communityId - 1],
                  graph.partition->partition.back());

        ///changing id of currentSingletonCommunity to last available id in current partition
        currentSingletonCommunity->communityId = int(graph.partition->partition.size());
    }


    static std::vector<Pair> getAggregatedNeighbours(CommunityStorage& community,
                                                     std::vector<pNode>& aggregatedGraph,
                                                     Graph& graph)    {
        auto aggregatedEdges = std::unordered_map<int, int>();

        for (auto& node : community) {
            for (auto& neighbourInfo : (*graph.adjacencyList)[node->vertexId - 1]) {
                const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);

                if (node->communityId != neighbour->communityId) {
                    if (aggregatedEdges.find(neighbour->communityId) == aggregatedEdges.end()) {
                        aggregatedEdges[neighbour->communityId] = neighbourInfo.edges;
                    } else {
                        aggregatedEdges[neighbour->communityId] += neighbourInfo.edges;
                    }
                }
            }
        }

        auto neighbours = std::vector<Pair>();
        for (const auto& element : aggregatedEdges) {
            auto pair = Pair(element.first, element.second);
            neighbours.push_back(pair);
        }

        return neighbours;
    }


    void aggregateNodes(Graph& graph) {
        std::vector<pNode> aggregatedNodes(graph.partition->partition.size());

        for (int i = 0; i < aggregatedNodes.size(); ++i) {
            auto aggregatedVertices = std::vector<pNode>();
            std::copy(graph.partition->partition[i]->community.begin(),
                      graph.partition->partition[i]->community.end(),
                      std::back_inserter(aggregatedVertices));

            auto vertex = std::make_shared<LouvainNode>(i + 1, 0, 0, 0);
            vertex->aggregatedVertices = std::move(aggregatedVertices);
            vertex->aggregationNumber = graph.partition->partition[i]->aggregationNumber;
            aggregatedNodes[i] = vertex;
        }

        auto aggregatedAdjacencyList = std::make_shared<std::vector<std::vector<Pair>>>(aggregatedNodes.size());
        for (int index = 0; index < aggregatedNodes.size(); ++index) {
            (*aggregatedAdjacencyList)[index] = getAggregatedNeighbours(graph.partition->partition[index]->community,
                                                                    aggregatedNodes, graph);
            aggregatedNodes[index]->loops = graph.partition->partition[index]->innerEdges;
            aggregatedNodes[index]->degree = graph.partition->partition[index]->communityDegree - 2 * aggregatedNodes[index]->loops;
        }

        graph.nodes = std::move(aggregatedNodes);
        graph.adjacencyList = std::move(aggregatedAdjacencyList);
    }


    std::vector<pCommunity> getSingletonPartition(Graph& graph) {
        std::vector<pCommunity> partition(graph.nodes.size());

        for (int index = 0; index < partition.size(); ++index) {
            graph.nodes[index]->communityId = index + 1;
            CommunityStorage community = {graph.nodes[index]};
            partition[index] = std::make_shared<Community>(community,
                                                           graph.nodes[index]->getFullDegree(),
                                                           graph.nodes[index]->loops);
            partition[index]->updateModularity(graph.edgesNumber);
            partition[index]->aggregationNumber =
                    std::dynamic_pointer_cast<LouvainNode>(graph.nodes[index])->aggregationNumber;
        }

        return partition;
    }


    static CommunityStorage flatCommunity(CommunityStorage& community, int communityId) {
        CommunityStorage flattenSet{};

        auto stack = std::vector<pNode>();
        std::copy(community.begin(), community.end(), std::back_inserter(stack));

        while (!stack.empty()) {
            auto louvainNode = std::dynamic_pointer_cast<LouvainNode>(stack.back());
            stack.pop_back();

            if (louvainNode->aggregatedVertices.empty()) {
                louvainNode->communityId = communityId;
                louvainNode->aggregationNumber = 1;
                flattenSet.insert(louvainNode);
                continue;
            }

            for (const auto& aggregatedVertex : louvainNode->aggregatedVertices) {
                stack.push_back(aggregatedVertex);
            }
        }

        return flattenSet;
    }

    std::vector<pCommunity> flat(const Partition& partition) {
        std::vector<pCommunity> flattenPartition(partition.partition.size());

        for (int index = 0; index < flattenPartition.size(); ++index) {
            auto flattenSet = flatCommunity(partition.partition[index]->community, index + 1);
            flattenPartition[index] = std::make_shared<Community>(flattenSet);
            flattenPartition[index]->aggregationNumber = flattenSet.size();
        }

        return flattenPartition;
    }
}