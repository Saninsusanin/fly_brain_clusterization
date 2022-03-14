#include "louvain_algorithm.h"
#include "../utils/utils.h"
#include "common.h"

#include <algorithm>
#include <unordered_map>
#include <memory>
#include <iostream>
#include <thread>

extern std::mt19937 generator;

namespace kv {
    static void insertNode(Partition& partition, pNode& node,
                           int maxCommunityId, double maxDelta, int edges,
                           std::unordered_map<int, int>& innerEdges) {
        node->communityId = maxCommunityId;
        partition.modularity += maxDelta;

        if (maxCommunityId > int(partition.partition.size())) {
            CommunityStorage newCommunity = {node};
            partition.partition.push_back(std::make_shared<Community>(newCommunity, node->getFullDegree(), node->loops));
            partition.partition.back()->updateModularity(edges);
            partition.partition.back()->aggregationNumber = std::dynamic_pointer_cast<LouvainNode>(node)->aggregationNumber;
        } else {
            partition.partition[maxCommunityId - 1]->insertNode(node, innerEdges[node->communityId], edges,
                                                                std::dynamic_pointer_cast<LouvainNode>(node)->aggregationNumber);
        }
    }

    static bool isInSingleton(Graph& graph, pNode& node) {
        if (graph.partition->partition[node->communityId - 1]->community.size() == 1) {
            swapCommunities(graph, node);
            return true;
        }

        return false;
    }


    ///deterministic move
    void deterministicMoveNode(pNode& node, Graph& graph) {
        ///check is currentNode in singleton community
        bool isAlreadyInSingleton = isInSingleton(graph, node);

        int maxCommunityId = node->communityId;
        auto innerEdges = node->getInnerEdgesFromNode(graph);
        double initialDelta = getInitialDelta(graph.partition->partition, node, isAlreadyInSingleton,
                                              graph.edgesNumber, innerEdges);
        double currentDelta = initialDelta;
        double maxDelta = 0.;

        ///iterating over all node's neighbours
        for (const auto& neighbourInfo : (*graph.adjacencyList)[node->vertexId - 1]) {
            const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
            currentDelta += getCurrentDelta(graph.partition->partition[neighbour->communityId - 1],
                                            node, graph.edgesNumber, innerEdges);
            ///remember current community id and value of metric if it is bigger that the previous one
            if (maxDelta < currentDelta) {
                maxDelta = currentDelta;
                maxCommunityId = neighbour->communityId;
            }

            currentDelta = initialDelta;
        }

        ///move to singleton community
        if (!isAlreadyInSingleton) {
            currentDelta += getCurrentSingletonDelta(node, graph.edgesNumber);

            if (maxDelta < currentDelta) {
                maxCommunityId = int(graph.partition->partition.size() + 1);
                maxDelta = currentDelta;
            }
        }

        ///move current node to community with max value of metric
        insertNode(*graph.partition, node, maxCommunityId, maxDelta, graph.edgesNumber, innerEdges);
    }


    ///stochastic move
    static void stochasticMoveNode(pNode& node, Graph& graph) {
        ///check is currentNode in singleton community
        bool isAlreadyInSingleton = isInSingleton(graph, node);

        int fitCommunityId = node->communityId;
        auto innerEdges = node->getInnerEdgesFromNode(graph);
        double initialDelta = getInitialDelta(graph.partition->partition, node, isAlreadyInSingleton,
                                              graph.edgesNumber, innerEdges);
        double currentDelta = initialDelta;
        double fitDelta = 0.;
        bool wasChanged = false;

        std::shuffle((*graph.adjacencyList)[node->vertexId - 1].begin(),
                     (*graph.adjacencyList)[node->vertexId - 1].end(),
                     generator);

        ///check all other communities
        for (const auto& neighbourInfo : (*graph.adjacencyList)[node->vertexId - 1]) {
            const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
            currentDelta += getCurrentDelta(graph.partition->partition[neighbour->communityId - 1],
                                            node, graph.edgesNumber, innerEdges);

            if (currentDelta > 0.) {
                fitDelta = currentDelta;
                fitCommunityId = neighbour->communityId;
                wasChanged = true;
                break;
            }
            currentDelta = initialDelta;
        }

        ///move to singleton community
        if (!isAlreadyInSingleton && !wasChanged) {
            currentDelta += getCurrentSingletonDelta(node, graph.edgesNumber);

            if (currentDelta > 0.) {
                fitCommunityId = int(graph.partition->partition.size() + 1);
                fitDelta = currentDelta;
            }
        }

        ///move current node to community with metric increasing
        insertNode(*graph.partition, node, fitCommunityId, fitDelta, graph.edgesNumber, innerEdges);
    }


    static void moveNodes(Graph& graph, std::function<void(pNode&, Graph&)>& moveNode) {
        double prevModularity;
        auto nodes(graph.nodes);

        std::shuffle(nodes.begin(), nodes.end(), generator);

        do {
            prevModularity = graph.partition->modularity;

            for (auto &&node : nodes) {
                moveNode(node, graph);
            }

        } while (prevModularity < graph.partition->modularity);
    }


    std::vector<pCommunity> clonePartitionData(pPartition& sourcePartition, Graph& targetGraph) {
        std::vector<pCommunity> partition(sourcePartition->partition.size());

        for (auto& currentPartition : sourcePartition->partition) {
            CommunityStorage community;

            for (auto& currentNode : currentPartition->community) {
                auto& targetNode = targetGraph.getNodeById(currentNode->vertexId);
                targetNode->communityId = currentNode->communityId;
                targetNode->loops = currentNode->loops;
                targetNode->degree = currentNode->degree;
                community.insert(targetNode);
            }

            partition[currentPartition->getCommunityId() - 1] =
                    std::make_shared<Community>(community,
                            currentPartition->communityDegree, currentPartition->innerEdges,
                            currentPartition->modularity, currentPartition->aggregationNumber);
        }

        return partition;
    }


    void Louvain(Graph& graph, bool notStochastic, bool notInitialized) {
        bool done;
        std::function<void(pNode&, Graph&)> moveNode = notStochastic ? deterministicMoveNode : stochasticMoveNode;

        std::shared_ptr<Partition> partition;
        Graph tmp;

        if (notInitialized) {
            tmp = graph;
            partition = std::make_shared<Partition>();
            partition->partition = getSingletonPartition(tmp);
            partition->setModularity();

            tmp.partition = partition;
        } else {
            tmp = *graph.copy();
            tmp.partition = std::make_shared<Partition>(graph.partition->modularity,
                    graph.partition->regularization,
                    clonePartitionData(graph.partition, tmp));
            tmp.partition->setModularity();

            partition = tmp.partition;
        }

        do {
            moveNodes(tmp, moveNode);
            done = (tmp.partition->partition.size() == tmp.nodes.size());

            if (!done) {
                aggregateNodes(tmp);
                partition->partition = getSingletonPartition(tmp);
                partition->setModularity();
            }

        }while(!done);

        partition->partition = flat(*partition);
        graph.partition = partition;
    }
}