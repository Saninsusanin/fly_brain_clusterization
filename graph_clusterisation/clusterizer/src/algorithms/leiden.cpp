#include <vector>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/random_device.hpp>

#include "leiden_algorithm.h"
#include "louvain_algorithm.h"
#include "common.h"


namespace kv {
    static void updateQueue(std::queue<pNode>& queue,
                            std::vector<int>& inQueue, const pNode& node,
                            const Graph& graph) {
        for (const auto& neighborInfo : (*graph.adjacencyList)[node->vertexId - 1]) {
            const auto& neighbor = graph.getNodeById(neighborInfo.vertexId);
            if (node->communityId != neighbor->communityId && inQueue[neighbor->vertexId] != 1) {
                queue.push(neighbor);
                inQueue[neighbor->vertexId] = 1;
            }
        }
    }

    static std::queue<pNode> fillQueue(std::vector<pNode>& nodes) {
        auto queue = std::queue<pNode>();

        for (const auto& vertex : nodes) {
            queue.push(vertex);
        }

        return queue;
    }

    void moveNodesFast(Graph& graph) {
        auto queue = fillQueue(graph.nodes);
        std::vector<int> inQueue(graph.nodes.size() + 1, 1);

        while (!queue.empty()) {
            auto node = queue.front();
            queue.pop();

            int prevCommunityId = node->communityId;
            inQueue[node->vertexId] = 0;

            deterministicMoveNode(node, graph);

            if (prevCommunityId != node->communityId) {
                updateQueue(queue, inQueue, node, graph);
            }
        }
    }

    static bool isInCommunity(const pNode& node, const Community& community) {
        return community.community.find(node) != community.community.end();
    }

///maybe to store number of common edges in Node structure
///@return amount of edges which are connecting
///@param node with current community
    static int getCommonEdges(const pNode& node, const Graph& graph, const Community& subset) {
        int commonEdges = 0;

        for (const auto& currentNeighbour : graph.getNeighbours(node->vertexId)) {
            commonEdges += isInCommunity(graph.getNodeById(currentNeighbour.vertexId), subset) ?
                           currentNeighbour.edges : 0;
        }

        return commonEdges;
    }

///@param graph - array of nodes in refined partition
///@param subset - community from partition given by MoveNodesFast function
///@return set of gamma-connected nodes
    static std::vector<pNode> fillInitialVertices(const Graph& graph, const Community& subset) {
        std::vector<pNode> initialVertices;

        for (const auto& currentNode : subset.community) {

            ///check if current node is gamma-connected with subset(in this case gamma is 1)
            if (getCommonEdges(currentNode, graph, subset) >=
                std::dynamic_pointer_cast<LouvainNode>(currentNode)->aggregationNumber *
                (subset.aggregationNumber - std::dynamic_pointer_cast<LouvainNode>(currentNode)->aggregationNumber)) {
                initialVertices.push_back(graph.nodes[currentNode->vertexId - 1]);
            }

        }

        return initialVertices;
    }

///get subcommunities of current subset
    static std::unordered_map<int, int> getSubcommunities(const Partition& partition, const Community& subset) {
        std::unordered_map<int, int> subCommunities;
        int numberOfVisited = 0;

        ///iterate over all communities in partition
        for (auto&& currentCommunity = partition.partition.begin();
             currentCommunity != partition.partition.end() && numberOfVisited < subset.community.size();
             ++currentCommunity) {
            auto subsetIter = subset.community.begin();
            auto currentCommunityIter = (*currentCommunity)->community.begin();

            while (currentCommunityIter != (*currentCommunity)->community.end() &&
                   subsetIter != subset.community.end()) {

                if ((*currentCommunityIter)->vertexId < (*subsetIter)->vertexId) {
                    break;
                } else if ((*currentCommunityIter)->vertexId > (*subsetIter)->vertexId) {
                    ++subsetIter;
                } else {
                    ++currentCommunityIter, ++subsetIter, ++numberOfVisited;

                    if (numberOfVisited == subset.community.size()) {
                        break;
                    }
                }

            }

            if (currentCommunityIter == (*currentCommunity)->community.end()) {
                subCommunities[(*(*currentCommunity)->community.begin())->communityId] =
                        (*(*currentCommunity)->community.begin())->communityId;
            }

        }

        return subCommunities;
    }

///getting common edges(in notation of Leiden paper - E(C, S - C))
    static int getCommonEdges(const Community& currentCommunity,
                              const Graph& graph,
                              const Community& subset) {
        int commonEdges = 0;

        for (const auto& currentNode : currentCommunity.community) {
            commonEdges += getCommonEdges(currentNode, graph, subset);
        }

        return commonEdges;
    }

    static bool isInSubcommunities(int communityId, std::unordered_map<int, int>& subcommunities) {
        return subcommunities.find(communityId) != subcommunities.end();
    }

///in notation of Leiden paper - E(C, S - C) >= ||C||*(||S||-||C||)
    static std::vector<int> getGammaConnectedCommunities(const pNode& node, const Graph& graph,
                                                         std::unordered_map<int, int>& subcommunities, const Community& subset) {
        std::vector<int> gammaConnectedCommunities;

        for (const auto& currentNeighbourInfo : graph.getNeighbours(node->vertexId)) {
            const auto& currentNeighbour = graph.getNodeById(currentNeighbourInfo.vertexId);
            auto currentCommunity = graph.partition->partition[currentNeighbour->communityId - 1];

            if (isInSubcommunities(currentNeighbour->communityId, subcommunities) &&
                getCommonEdges(*currentCommunity, graph, subset) >= currentCommunity->aggregationNumber *
                                                                    (subset.aggregationNumber - currentCommunity->aggregationNumber)) {
                gammaConnectedCommunities.push_back(currentNeighbour->communityId);
            }

        }

        return gammaConnectedCommunities;
    }

///@param node - node from partition given by refinedPartition
///@return is in singleton community in refined partition
    static bool isInSingleton(const pNode& node, const Graph& graph) {
        return graph.partition->partition[node->communityId - 1]->community.size() == 1;
    }

    using Deltas = std::vector<std::pair<int, double>>;

    static Deltas getAllDeltas(Partition& partition, std::vector<int>& gammaConnectedCommunities,
                               std::unordered_map<int, int>& innerEdges, const pNode& node, int numberOfEdges) {
        Deltas deltas(gammaConnectedCommunities.size());
        auto initialDelta = getInitialDelta(partition.partition, node,
                                            true, numberOfEdges, innerEdges);

        for (int currentCommunityIndex = 0;
             currentCommunityIndex < gammaConnectedCommunities.size();
             ++currentCommunityIndex) {
            auto currentDelta = getCurrentDelta(partition.partition[gammaConnectedCommunities[currentCommunityIndex]],
                                                node, numberOfEdges, innerEdges);
            deltas[currentCommunityIndex] = std::pair<int, double>(gammaConnectedCommunities[currentCommunityIndex],
                                                                   initialDelta + currentDelta);
        }

        return deltas;
    }

    static std::vector<double> getNonStandartized(Deltas& deltas, bool& isAllZeros) {
        std::vector<double> nonStandardized(deltas.size());

        for (int currentDeltaIndex = 0; currentDeltaIndex < deltas.size(); ++currentDeltaIndex) {

            if (deltas[currentDeltaIndex].second > 0 + 1e-14) {
                isAllZeros = false;
                nonStandardized[currentDeltaIndex] =
                        exp(deltas[currentDeltaIndex].second / LeidenConfig::TETHA);
            }

        }

        return nonStandardized;
    }

    static int getRandomCommunityIndex(std::vector<double>& probabilities) {
        boost::random::random_device rand_dev;
        boost::mt19937 gen(rand_dev());
        boost::random::discrete_distribution<> dist(probabilities);

        return dist(gen);
    }

    static void updateSubcommunities(std::unordered_map<int, int>& subcommunities,
                                     const Graph& graph, const pNode& node) {
        if (!isInSubcommunities(graph.partition->partition.size(), subcommunities)) {
            subcommunities.erase(node->communityId);
            subcommunities[graph.partition->partition.size()] = graph.partition->partition.size();
        }
    }

    void mergeNodesSubset(Graph& graph, const Community& subset) {
        int destinationCommunityIndex;

        auto initialVertices = fillInitialVertices(graph, subset);
        auto subcommunities = getSubcommunities(*(graph.partition), subset);

        for (auto& node : initialVertices) {

            if (isInSingleton(node, graph)) {
                auto innerEdges = node->getInnerEdgesFromNode(graph);
                updateSubcommunities(subcommunities, graph, node);
                swapCommunities(graph, node);
                auto gammaConnectedCommunities = getGammaConnectedCommunities(node, graph,
                                                                              subcommunities, subset);
                auto deltas = getAllDeltas(*(graph.partition), gammaConnectedCommunities,
                                           innerEdges, node, graph.edgesNumber);
                bool isAllZeros = true;
                auto nonStandardized = getNonStandartized(deltas, isAllZeros);

                ///if isAllZeros is true this is mean that current vertex must stay in singleton community
                if (isAllZeros) {
                    node->communityId = (int)graph.partition->partition.size() + 1;
                    graph.partition->partition.push_back(std::make_shared<Community>(
                            Community(CommunityStorage {node},
                                      node->getFullDegree(), node->loops, 0,
                                      std::dynamic_pointer_cast<LouvainNode>(node)->aggregationNumber)));
                    graph.partition->partition.back()->updateModularity(graph.edgesNumber);
                } else {
                    subcommunities.erase(node->communityId);
                    destinationCommunityIndex = getRandomCommunityIndex(nonStandardized);
                    node->communityId = deltas[destinationCommunityIndex].first;
                    graph.partition->partition[deltas[destinationCommunityIndex].first - 1]->community.insert(node);
                    graph.partition->partition[deltas[destinationCommunityIndex].first - 1]->innerEdges +=
                            innerEdges[deltas[destinationCommunityIndex].first];
                    graph.partition->partition[deltas[destinationCommunityIndex].first - 1]->aggregationNumber +=
                            std::dynamic_pointer_cast<LouvainNode>(node)->aggregationNumber;
                    graph.partition->partition[deltas[destinationCommunityIndex].first - 1]->communityDegree +=
                            node->getFullDegree();
                    graph.partition->partition[deltas[destinationCommunityIndex].first - 1]->updateModularity(graph.edgesNumber);
                    graph.partition->modularity += deltas[destinationCommunityIndex].second;
                }
            }
        }
    }

    void refinePartition(Graph& refineGraph, const Graph& moveGraph) {
        auto singletonPartition = getSingletonPartition(refineGraph);
        auto refinePartition = std::make_shared<Partition>(0, 0, singletonPartition);
        refinePartition->setMetricValue(refineGraph.nodes.size(), refineGraph.edgesNumber);
        refineGraph.partition = refinePartition;

        for (auto& currentCommunity : moveGraph.partition->partition) {
            mergeNodesSubset(refineGraph, *currentCommunity);
            //std::cout << (*currentCommunity->community.begin())->communityId << std::endl;
        }

        //std::cout << "refined";
    }

    void updateMovePartition(const Graph& refineGraph, Graph& moveGraph) {
        auto newGraph = std::make_shared<Graph>(std::vector<pNode>(refineGraph.nodes.size()),
                                                refineGraph.adjacencyList,
                                                refineGraph.edgesNumber,
                                                std::make_shared<Partition>(0, 0,
                                                                            std::vector<pCommunity>(moveGraph.partition->partition.size(), nullptr)));

        ///iterate over all aggregated nodes in refined graph
        for (int currentNodeId = 1; currentNodeId <= refineGraph.nodes.size(); ++currentNodeId) {
            auto LouvainCurrentNode = std::dynamic_pointer_cast<LouvainNode>(refineGraph.getNodeById(currentNodeId));
            ///initialize new aggregated vertex in moved graph
            auto aggregatedNode = std::make_shared<LouvainNode>(currentNodeId,
                                                                moveGraph.getNodeById(LouvainCurrentNode->aggregatedVertices[0]->vertexId)->communityId,
                                                                LouvainCurrentNode->degree, LouvainCurrentNode->loops,
                                                                std::vector<pNode> (LouvainCurrentNode->aggregatedVertices.size()),
                                                                LouvainCurrentNode->aggregationNumber);
            newGraph->nodes[currentNodeId - 1] = aggregatedNode;

            ///remember all vertices from which current aggregated node was aggregated
            for (int aggregatedIndex = 0;
                 aggregatedIndex < LouvainCurrentNode->aggregatedVertices.size();
                 ++aggregatedIndex) {
                aggregatedNode->aggregatedVertices[aggregatedIndex] =
                        moveGraph.getNodeById(LouvainCurrentNode->aggregatedVertices[aggregatedIndex]->vertexId);
            }

            ///add current aggregated node to the suitable community
            if (newGraph->partition->partition[aggregatedNode->communityId - 1] == nullptr) {
                auto tmpCommunity = moveGraph.partition->partition
                [moveGraph.getNodeById(LouvainCurrentNode->aggregatedVertices[0]->vertexId)->communityId - 1];
                newGraph->partition->partition[aggregatedNode->communityId - 1] =
                        std::make_shared<Community>(CommunityStorage{aggregatedNode},
                                                    tmpCommunity->communityDegree, tmpCommunity->innerEdges,
                                                    tmpCommunity->modularity, tmpCommunity->aggregationNumber);
            } else {
                newGraph->partition->partition[aggregatedNode->communityId - 1]->community.insert(aggregatedNode);
            }

        }

        newGraph->partition->modularity = moveGraph.partition->modularity;
        moveGraph.nodes = std::move(newGraph->nodes);
        moveGraph.partition = newGraph->partition;
        moveGraph.adjacencyList = newGraph->adjacencyList;
    }

    static std::vector<pNode> getNodesFromPartition(const pPartition& partition, int numberOfVertices) {
        std::vector<pNode> nodes(numberOfVertices);

        for (const auto& currentCommunity : partition->partition) {

            for (const auto& currentNode : currentCommunity->community) {
                nodes[currentNode->vertexId - 1] = currentNode;
            }

        }

        return nodes;
    }

    void Leiden(Graph& graph) {
        bool done = false;
        auto copyMoveGraph = graph.copy();
        auto partition = std::make_shared<Partition>(0, 0, getSingletonPartition(*copyMoveGraph));
        partition->setModularity();
        copyMoveGraph->partition = partition;
        auto copyRefineGraph = graph.copy();
        int numberOfIterations = 0;

        while(!done) {
            moveNodesFast(*copyMoveGraph);
            refinePartition(*copyRefineGraph, *copyMoveGraph);
            aggregateNodes(*copyRefineGraph);
            copyRefineGraph->partition->partition = getSingletonPartition(*copyRefineGraph);
            copyRefineGraph->partition->setModularity();
            updateMovePartition(*copyRefineGraph, *copyMoveGraph);
            done = copyRefineGraph->partition->partition.size() == copyRefineGraph->nodes.size();
            ++numberOfIterations;
        }

        graph.partition = std::make_shared<Partition>(0, 0, std::vector<pCommunity>());
        graph.partition->partition = flat(*(copyMoveGraph->partition));
        graph.partition->modularity = copyMoveGraph->partition->modularity;
        graph.nodes = getNodesFromPartition(graph.partition, graph.nodes.size());
        //std::cout << numberOfIterations << std::endl;
    }
}