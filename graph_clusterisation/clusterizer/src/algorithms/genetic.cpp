#include "genetic_algorithm.h"
#include "../utils/utils.h"

#include <vector>
#include <algorithm>
#include <memory>
#include <future>
#include <array>
#include <cmath>
#include <functional>
#include <random>
#include <iostream>
#include <iomanip>
#include <boost/functional/hash.hpp>


extern std::mt19937 generator;


namespace kv {
    static void updateGraph(std::vector<pNode>& nodes) {
        for (auto& node : nodes) {
            node->communityId = 0;
        }

    }

    static void connectedTraverse(std::queue<pNode>& queue,
                                  std::vector<pCommunity>& partitionData,
                                  Graph& graph) {
        while (!queue.empty()) {
            auto currentNode = queue.front();

            queue.pop();

            auto neighbours = graph.getNeighbours(currentNode->vertexId);
            for (auto& neighbourInfo : neighbours) {
                const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);

                if (neighbour->communityId == 0) {
                    neighbour->communityId = currentNode->communityId;
                    queue.push(neighbour);
                    partitionData[neighbour->communityId - 1]->community.insert(neighbour);
                }
            }
        }
    }

    static void disconnectedTraverse(std::queue<pNode>& queue,
                                     std::vector<pCommunity>& partitionData,
                                     std::vector<pNode>& nodes, int& numberOfCommunities,
                                     Graph& graph) {
        auto&& vertex = nodes.begin();

        while (vertex != nodes.end()) {
            ///find new community center
            for (;vertex != nodes.end(); ++vertex) {

                if ((*vertex)->communityId == 0) {
                    (*vertex)->communityId = ++numberOfCommunities;
                    CommunityStorage currentCommunity = {*vertex};
                    partitionData.push_back(std::make_shared<Community>(currentCommunity));
                    queue.push(*vertex);

                    break;
                }

            }

            connectedTraverse(queue, partitionData, graph);
        }
    }

    void kNearestNeighboursPartition(Graph& graph) {
        auto nodes(graph.nodes);
        updateGraph(nodes);

        ///sort graph by id of vertex
        std::sort(nodes.begin(), nodes.end(),
                  [](const pNode& left, const pNode& right)
                  { return (left->degree > right->degree); });

        ///select number of vertices from which center of community would be selected
        int numberOfShuffledElements = nodes[0]->degree;

        ///shuffle of first numberOfShuffledElements
        std::shuffle(nodes.begin(), nodes.begin() + numberOfShuffledElements - 1, generator);

        ///getting random number of communities
        int numberOfCommunities = int(Utils::generateRandomNumber(kNearestConfig::MIN_NUMBER_OF_COMMUNITIES,
                                                                  kNearestConfig::MAX_NUMBER_OF_COMMUNITIES));
        std::vector<pCommunity> partitionData;
        std::queue<pNode> queue;

        ///generate first numberOfCommunities communities
        for (int index = 0; index < numberOfCommunities; ++index) {
            nodes[index]->communityId = index + 1;
            CommunityStorage currentCommunity = {nodes[index]};
            partitionData.push_back(std::make_shared<Community>(currentCommunity));
            queue.push(nodes[index]);
        }

        ///traverse all nodes in current linked component
        connectedTraverse(queue, partitionData, graph);

        ///traverse nodes which is in another connectivity component
        disconnectedTraverse(queue, partitionData, nodes, numberOfCommunities, graph);

        auto partition = std::make_shared<Partition>();
        partition->partition = std::move(partitionData);
        partition->updatePartition(graph);
        partition->setMetricValue(graph.nodes.size(), graph.edgesNumber);

        graph.partition = partition;
    }


    void randomInitialization(Graph& graph, int numberOfClusters) {
        double step = graph.nodes.size() / double(numberOfClusters);
        std::shuffle(graph.nodes.begin(), graph.nodes.end(), generator);

        std::vector<pCommunity> partitionData;

        for (int clusterID = 0; clusterID < numberOfClusters; ++clusterID) {

            for (int index = int(clusterID * step); index < int((clusterID + 1) * step); ++index) {
                graph.nodes[index]->communityId = clusterID + 1;
            }

            CommunityStorage currentCommunity(graph.nodes.begin() + int(clusterID * step),
                    graph.nodes.begin() + int((clusterID + 1) * step));
            partitionData.push_back(std::make_shared<Community>(currentCommunity));
        }

        auto partition = std::make_shared<Partition>();
        partition->partition = std::move(partitionData);
        partition->updatePartition(graph);
        partition->setMetricValue(graph.nodes.size(), graph.edgesNumber);
        std::sort(graph.nodes.begin(), graph.nodes.end(),
                [](pNode& left, pNode& right) {return left->vertexId < right->vertexId;});

        graph.partition = partition;
    }


    static CommunityStorage intersectSets(CommunityStorage& firstSet,
                                          CommunityStorage& secondSet,
                                          Graph& thirdGraph,
                                          int currentCommunityIndex) {
        CommunityStorage intersection{};
        auto firstIter = firstSet.begin();
        auto secondIter = secondSet.begin();

        while (firstIter != firstSet.end() && secondIter != secondSet.end()) {
            if ((*firstIter)->vertexId == (*secondIter)->vertexId) {
                const auto& newNode = thirdGraph.getNodeById((*firstIter)->vertexId);
                newNode->communityId = currentCommunityIndex;
                intersection.insert(newNode);
                firstIter++;
                secondIter++;
            } else if ((*firstIter)->vertexId > (*secondIter)->vertexId) {
                secondIter++;
            } else {
                firstIter++;
            }
        }

        return intersection;
    }


    static void setCommunityId(CommunityStorage& community, int communityId) {
        for (auto& node : community) {
            node->communityId = communityId;
        }
    }


    pGraph cross(const Graph& firstGraph, const Graph& secondGraph) noexcept {
        auto thirdGraph = firstGraph.copy();

        std::vector<pCommunity> newPartitionData{};
        int currentCommunityIndex = 1;

        for (const auto& firstIter : firstGraph.partition->partition) {
            for (const auto& secondIter : secondGraph.partition->partition)
            {
                auto intersection = intersectSets(firstIter->community, secondIter->community,
                                                  *thirdGraph, currentCommunityIndex);

                if (!intersection.empty()) {
                    auto intersectionCommunity = std::make_shared<Community>(intersection);
                    newPartitionData.push_back(intersectionCommunity);
                    currentCommunityIndex++;
                }
            }
        }

        auto partition = std::make_shared<Partition>();
        partition->partition = newPartitionData;
        partition->updatePartition(*thirdGraph);
        partition->setMetricValue(thirdGraph->getVerticesNumber(), thirdGraph->edgesNumber);

        thirdGraph->partition = std::move(partition);

        return thirdGraph;
    }


/// Mutation is performed by merging some communities
    void mergingMutate(Graph& graph, GeneticConfig& config) noexcept {
        if (graph.partition->partition.size() == 1) {
            return;
        }

        std::shuffle(graph.partition->partition.begin(),
                     graph.partition->partition.end(),
                     generator);

        /// maximal number of merging communities should not be less than minimal
        int maxMergingCommunities = std::max(int(config.MAX_MERGING_PART *
                                                 graph.partition->partition.size()),
                                             int(config.MIN_MERGING_COMMUNITIES));

        int mergingNumber = int(Utils::generateRandomNumber(config.MIN_MERGING_COMMUNITIES,
                                                            maxMergingCommunities));

        /// merging <mergingNumber> last communities
        CommunityStorage mergedCommunities{};

        for (int iteration = 0; iteration < mergingNumber; ++iteration) {
            auto& lastCommunity = graph.partition->partition.back();
            for (const auto& node : lastCommunity->community) {
                mergedCommunities.insert(node);
            }

            graph.partition->partition.pop_back();
        }


        graph.partition->partition.push_back(std::make_shared<Community>(mergedCommunities));

        for (int communityId = 1; communityId <= graph.partition->partition.size(); communityId++) {
            setCommunityId(graph.partition->partition[communityId - 1]->community, communityId);
        }

        graph.partition->partition.back()->updateCommunityDegree();
        graph.partition->partition.back()->updateInnerEdges(graph);
        graph.partition->setMetricValue(graph.nodes.size(), graph.edgesNumber);

    }


    static CommunityStorage::iterator getMinDegreeNode(CommunityStorage& community) {
        auto minDegreeNode = community.begin();
        for (auto iter = community.begin(); iter != community.end(); iter++) {
            if ((*iter)->degree < (*minDegreeNode)->degree) {
                minDegreeNode = iter;
            }
        }
        return minDegreeNode;
    }


/// Mutation is performed by extracting vertices from
/// the least dense communities to a singleton communities
    void extractingMutate(Graph& graph, GeneticConfig& config) noexcept {
        /// sort communities by density value (in increasing order)
        std::sort(graph.partition->partition.begin(),
                  graph.partition->partition.end(),
                  [](const pCommunity& a, const pCommunity& b)
                  {return a->getCommunityDensity() < b->getCommunityDensity();});

        int maxCommunities = std::max(int(config.MIN_EXTRACTING_COMMUNITIES),
                                      int(config.MAX_EXTRACTING_COMMUNITIES_PART * graph.partition->partition.size()));

        int numberOfCommunities = int(Utils::generateRandomNumber(config.MIN_EXTRACTING_COMMUNITIES,
                                                              maxCommunities));

        for (int communityIndex = 0; communityIndex < numberOfCommunities; ++communityIndex) {
            if (graph.partition->partition[communityIndex]->community.size() == 1) {
                continue;
            }

            int maxNodes = std::max(int(config.MIN_EXTRACTING_NODES),
                    int(config.MAX_EXTRACTING_NODES_PART * graph.partition->partition[communityIndex]->community.size()));

            int numberOfVertices = int(Utils::generateRandomNumber(config.MIN_EXTRACTING_NODES, maxNodes));

            for (int vertexIndex = 0; vertexIndex < numberOfVertices; ++vertexIndex) {
                auto extractedNode = getMinDegreeNode(graph.partition->partition[communityIndex]->community);
                (*extractedNode)->communityId = int(graph.partition->partition.size()) + 1;

                CommunityStorage singletonCommunity = {*extractedNode};
                graph.partition->partition.push_back(std::make_shared<Community>(singletonCommunity));
                graph.partition->partition[communityIndex]->community.erase(extractedNode);
            }
        }

        graph.partition->updatePartition(graph);
        graph.partition->setMetricValue(graph.nodes.size(), graph.edgesNumber);
    }

    /// plug for Girvan-Newman algorithm
    std::vector<std::vector<int>> randomSplit(std::vector<pNode>& community) noexcept {
        std::vector<std::vector<int>> newCommunities;
        std::shuffle(community.begin(), community.end(), generator);

        newCommunities.emplace_back();
        newCommunities.emplace_back();

        int i = 0;

        for (auto& node : community) {
            if (i < community.size() / 2) {
                newCommunities[0].push_back(node->vertexId);
                ++i;
            } else {
                newCommunities[1].push_back(node->vertexId);
            }
        }

        return newCommunities;
    }


    /// current version is not working with loops in a graph(maybe)
    std::vector<std::vector<int>> GirvanNewman(Graph& graph, std::vector<pNode>& community) noexcept {
        std::vector<std::vector<int>> newCommunities;
        std::unordered_map<std::pair<int, int>, std::shared_ptr<EdgeInfo>, boost::hash<std::pair<int, int>>> edgesInfo;
        std::queue<pNode> queue;
        std::unordered_map<int, int> indexHolder;
        std::vector<NodeInfo> nodesInfo(community.size());
        std::vector<int> indices(community.size());
        bool done;
        size_t sum;

        /// association of the indices in the isUsed and the numberOfPaths vectors with vertexID
        for (int indexInCommunity = 0; indexInCommunity < community.size(); ++indexInCommunity) {
            indexHolder.insert({community[indexInCommunity]->vertexId, indexInCommunity});
        }

        /// initialize edgesInfo
        for (auto& node : community) {
            auto neighbours = graph.getNeighbours(node->vertexId);

            for (auto& neighbourInfo : neighbours) {
                const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);

                if (neighbour->communityId == node->communityId) {
                    edgesInfo.insert({{node->vertexId < neighbourInfo.vertexId ?
                                              node->vertexId : neighbourInfo.vertexId,
                                              node->vertexId < neighbourInfo.vertexId ?
                                              neighbourInfo.vertexId : node->vertexId},
                                      std::make_shared<EdgeInfo>()});
                }
            }
        }

        /// delete edges with the highest edge betweenness until a graph becomes disconnected
        do {
            sum = 0;

            /// reset edgesBetweenness
            for (auto& currentInfo : edgesInfo) {
                currentInfo.second->edgeBetweenness = 0;
            }

            for (int indexInCommunity = 0; indexInCommunity < community.size(); ++indexInCommunity) {
                /// set statuses of each node
                /// initialize number of current paths
                /// initialize level of each node in current traverse
                for (int i = 0; i < community.size(); ++i) {
                    nodesInfo[i] = {false,0,
                                    0, 1., community[i]->vertexId};
                }

                /// preparation of the current center of the traverse
                nodesInfo[indexInCommunity] =
                        {true,1, 0,
                         1.,community[indexInCommunity]->vertexId};
                queue.push(community[indexInCommunity]);

                /// bfs
                while (!queue.empty()) {
                    auto& currentNode = queue.front();
                    queue.pop();
                    auto neighbours = graph.getNeighbours(currentNode->vertexId);

                    for (auto& neighbourInfo : neighbours) {
                        const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
                        int index = indexHolder[neighbour->vertexId];

                        if (neighbour->communityId == community[0]->communityId) {
                            auto& edge = edgesInfo[{currentNode->vertexId < neighbourInfo.vertexId ?
                                                    currentNode->vertexId : neighbourInfo.vertexId,
                                                    currentNode->vertexId < neighbourInfo.vertexId ?
                                                    neighbourInfo.vertexId : currentNode->vertexId}];

                            if (!edge->isDeleted) {

                                if (!nodesInfo[index].isUsed) {
                                    nodesInfo[index].isUsed = true;
                                    nodesInfo[index].level = nodesInfo[indexHolder[currentNode->vertexId]].level + 1;
                                    queue.push(neighbour);
                                }

                                if (nodesInfo[index].level > nodesInfo[indexHolder[currentNode->vertexId]].level) {
                                    nodesInfo[index].numberOfPaths +=
                                            nodesInfo[indexHolder[currentNode->vertexId]].numberOfPaths;
                                }
                            }
                        }


                    }
                }
                /// the section below updates edgeBetweenness values
                /// fill index array
                int i = 0;

                for (auto& index : indices) {
                    index = i++;
                }

                std::sort(indices.begin(), indices.end(),
                        [nodesInfo](int& left, int& right) {
                    return nodesInfo[left].level > nodesInfo[right].level;
                });

                /// actual edge betweenness update
                for (auto index : indices) {
                    auto neighbours = graph.getNeighbours(nodesInfo[index].vertexId);

                    for (auto& neighbourInfo : neighbours) {
                        const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);

                        if (neighbour->communityId == community[0]->communityId) {
                            auto& edge = edgesInfo[{nodesInfo[index].vertexId < neighbourInfo.vertexId ?
                                                    nodesInfo[index].vertexId : neighbourInfo.vertexId,
                                                    nodesInfo[index].vertexId < neighbourInfo.vertexId ?
                                                    neighbourInfo.vertexId : nodesInfo[index].vertexId}];

                            if (!edge->isDeleted) {
                                const auto neighbourIndex = indexHolder[neighbour->vertexId];
                                double proportion = nodesInfo[neighbourIndex].numberOfPaths /
                                                    double(nodesInfo[index].numberOfPaths);

                                if (nodesInfo[index].level > nodesInfo[neighbourIndex].level) {
                                    nodesInfo[neighbourIndex].ebStorage += proportion;
                                    edge->edgeBetweenness = nodesInfo[index].ebStorage * proportion;
                                }
                            }
                        }
                    }
                }
            }

            /// delete edges with the highest edgeBetweenness value
            std::vector<std::shared_ptr<EdgeInfo>> tmp;
            tmp.reserve(edgesInfo.size());
            std::transform(edgesInfo.begin(), edgesInfo.end(), std::back_inserter(tmp),
                    [](std::pair<const std::pair<int, int>, std::shared_ptr<EdgeInfo>>& elem) { return elem.second; });
            std::sort(tmp.begin(), tmp.end(),
                    [](std::shared_ptr<EdgeInfo>& left, std::shared_ptr<EdgeInfo>& right) {
                return left->edgeBetweenness > right->edgeBetweenness;
            });

            for (int i = 0; i < tmp.size(); ++i) {
                if (std::fabs(tmp[i]->edgeBetweenness - tmp[0]->edgeBetweenness) < 1e-14) {
                    tmp[i]->isDeleted = true;
                } else {
                    break;
                }
            }

            for (const auto& status : nodesInfo) { sum += status.isUsed; }

            /// check whether graph is connected or not
            done = sum != newCommunities.size();
        }while (!done);

        /// iterate over first connected component
        std::vector<int> newCommunityIndices;

        for (int i = 0; i < community.size(); ++i) {
            if (nodesInfo[i].isUsed) {
                newCommunityIndices.push_back(nodesInfo[i].vertexId);
            }
        }

        newCommunities.push_back(newCommunityIndices);

        while (sum < nodesInfo.size()) {
            std::vector<int> tmpCommunityStorage;
            /// insert new node to the queue
            for (int i = 0; i < community.size(); ++i) {
                if (!nodesInfo[i].isUsed) {
                    nodesInfo[i].isUsed = true;
                    ++sum;
                    tmpCommunityStorage.push_back(nodesInfo[i].vertexId);
                    queue.push(community[i]);
                }
            }

            while (!queue.empty()) {
                auto& currentNode = queue.front();
                queue.pop();
                auto neighbours = graph.getNeighbours(currentNode->vertexId);

                for (auto& neighbourInfo : neighbours) {
                    const auto& neighbour = graph.getNodeById(neighbourInfo.vertexId);
                    int index = indexHolder[neighbour->vertexId];

                    if (neighbour->communityId == community[0]->communityId &&
                        !edgesInfo[{currentNode->vertexId < neighbourInfo.vertexId ?
                                    currentNode->vertexId : neighbourInfo.vertexId,
                                    currentNode->vertexId < neighbourInfo.vertexId ?
                                    neighbourInfo.vertexId : currentNode->vertexId}]->isDeleted &&
                                    !nodesInfo[index].isUsed) {
                        nodesInfo[index].isUsed = true;
                        ++sum;
                        tmpCommunityStorage.push_back(nodesInfo[index].vertexId);
                        queue.push(neighbour);
                    }
                }
            }

            newCommunities.push_back(tmpCommunityStorage);
        }

        return newCommunities;
    }


    void separatingMutate(Graph& graph, GeneticConfig& config,
            std::function<std::vector<std::vector<int>>(std::vector<kv::pNode>&)>& engine) noexcept {
        std::shuffle(graph.partition->partition.begin(), graph.partition->partition.end(), generator);

        int maxCommunities = std::max(int(config.MIN_SEPARATING_COMMUNITIES),
                                      int(config.MAX_SEPARATING_COMMUNITIES_PART * graph.partition->partition.size()));
        int numberOfCommunities = int(Utils::generateRandomNumber(config.MIN_SEPARATING_COMMUNITIES,
                                                                  maxCommunities));
        int currentMaxCommunityID = graph.partition->partition.size();

        for (int communityIndex = 0; communityIndex < numberOfCommunities; ++communityIndex) {
            auto community = graph.partition->partition[communityIndex]->community;

            if (community.size() == 1) {
                continue;
            }

            std::vector<pNode> communityVec;
            std::copy(community.begin(), community.end(), std::back_inserter(communityVec));
            std::vector<std::vector<int>> newCommunities = engine(communityVec);

            for (int newCommunityID = 0; newCommunityID < newCommunities.size(); ++newCommunityID) {
                std::vector<pNode> newCommunityVec;

                for (auto vertexID : newCommunities[newCommunityID]) {
                    const auto& node = graph.getNodeById(vertexID);
                    node->communityId = newCommunityID == 0 ? node->communityId : currentMaxCommunityID + newCommunityID;
                    newCommunityVec.push_back(node);
                }

                CommunityStorage newCommunity(
                        std::make_move_iterator(newCommunityVec.begin()),
                        std::make_move_iterator(newCommunityVec.end()));

                if (newCommunityID == 0) {
                    graph.partition->partition[communityIndex] = std::make_shared<Community>(newCommunity);
                } else {
                    graph.partition->partition.push_back(std::make_shared<Community>(newCommunity));
                }

            }

            currentMaxCommunityID += int(newCommunities.size() - 1);
        }

        graph.partition->updatePartition(graph);
        graph.partition->setMetricValue(graph.nodes.size(), graph.edgesNumber);
    }


    using indices_t = std::vector<std::array<int, 3>>;
    using population_iter_t = std::vector<pGraph>::iterator;


    ///pairwise cross
    static void getNewPopulation(std::vector<pGraph>& population,
            indices_t::iterator begin,
            indices_t::iterator end) {
        for (; begin != end; ++begin){
            population[(*begin)[Indices::PLACE_ID]] = cross(
                    *population[(*begin)[Indices::FIRST_PARENT_ID]],
                    *population[(*begin)[Indices::SECOND_PARENT_ID]]);
        }
    }


    static void getFirstPopulation(Graph& graph,
            population_iter_t begin,
            population_iter_t end,
            std::function<void(Graph&)>& method){
        for (;begin != end; ++begin) {
            *begin = graph.copy();
            method(**begin);
        }
    }


    static void fillParentIndices(indices_t& indices, int populationSize) {
        for (int first = 0, index = 0; first < populationSize; ++first) {
            for (int second = first + 1; second < populationSize; ++second, ++index) {
                indices[index] = {first, second, populationSize + index};
            }
        }
    }


    static void generateFirstPopulation(Graph& graph, std::vector<std::future<void>>& futures,
                                          std::vector<pGraph>& population, int sizeOfPopulation,
                                          GeneticConfig& config, std::function<void(Graph&)>& method) {
        int denominator = std::min(config.NUMBER_OF_THREADS, sizeOfPopulation);
        double step = sizeOfPopulation / double(denominator);
        futures.reserve(denominator);
        futures.resize(futures.capacity());

        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID] = std::move(std::async(
                    std::launch::async,
                    ///worker
                    getFirstPopulation,
                    ///reference for current graph
                    std::ref(graph),
                    ///begin of the interest region for this thread
                    population.begin() + int(step * threadID),
                    ///end of the interest region for this thread
                    population.begin() + int(step * (threadID + 1)),
                    ///the first population generation method
                    std::ref(method)));
        }

        ///thread join
        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID].get();
        }
    }


    static void applyCross(std::vector<std::future<void>>& futures, indices_t& parentIndices,
                            std::vector<pGraph>& population, int denominator,
                            double step, int amountOfNew) {
        ///split work in several threads
        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID] = std::move(std::async(
                    std::launch::async,
                    ///worker
                    getNewPopulation,
                    ///reference to the population
                    std::ref(population),
                    ///iterator on the begin of the region to work with
                    parentIndices.begin() + int(threadID * step),
                    ///iterator on the end of current amount of indices to work with
                    parentIndices.begin() + int(step * (threadID + 1))));
        }

        ///thread join
        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID].get();
        }
    }


    static void greedyRefinement(population_iter_t begin, population_iter_t end) {
        for (;begin != end; ++begin) {
            Louvain(**begin, false, false);
        }
    }


    static void refinePopulation(std::vector<std::future<void>>& futures,
                                  std::vector<pGraph>& population, int sizeOfPopulation,
                                  GeneticConfig& config) {
        int amountOfNew = sizeOfPopulation * (sizeOfPopulation - 1) / 2;
        int denominator = std::min(config.NUMBER_OF_THREADS, amountOfNew);
        double step = amountOfNew / double(denominator);
        futures.reserve(denominator);
        futures.resize(futures.capacity());

        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID] = std::move(std::async(
                    std::launch::async,
                    ///worker
                    greedyRefinement,
                    ///begin of the interest region for this thread
                    population.begin() + sizeOfPopulation + int(step * threadID),
                    ///end of the interest region for this thread
                    population.begin() + sizeOfPopulation + int(step * (threadID + 1))));
        }

        ///thread join
        for (int threadID = 0; threadID < denominator; ++threadID) {
            futures[threadID].get();
        }
    }


    void geneticAlgorithm(Graph& graph, std::function<void(Graph&)>& method,
                          std::vector<std::function<void(kv::Graph&, GeneticConfig&)>>& mutations,
                          GeneticConfig& config) {
        int sizeOfPopulation = int(Utils::generateRandomNumber(config.MIN_SIZE_OF_POPULATION,
                                                               config.MAX_SIZE_OF_POPULATION));
        //BOOST_LOG_TRIVIAL(debug) << "size of the population: " << sizeOfPopulation;
        std::cout << "size of the population: " << sizeOfPopulation << '\n';
        int amountOfNew = sizeOfPopulation * (sizeOfPopulation - 1) / 2;
        int current_schedule = 0;
        double previousMetricValue;

        std::vector<pGraph> population;
        population.reserve(sizeOfPopulation + amountOfNew);
        population.resize(population.capacity());

        ///storage of futures
        std::vector<std::future<void>> futures;

        ///storage of data for each thread in cross
        indices_t parentIndices;
        parentIndices.reserve(amountOfNew);
        parentIndices.resize(parentIndices.capacity());
        fillParentIndices(parentIndices, sizeOfPopulation);

        ///generation of first population
        generateFirstPopulation(graph, futures, population, sizeOfPopulation, config, method);

        ///reset step for the parallelization of crosses
        int denominator = std::min(config.NUMBER_OF_THREADS, amountOfNew);
        double step = amountOfNew / double(denominator);
        futures.reserve(denominator);
        futures.resize(futures.capacity());

        for (int iteration = 0; iteration < config.MAX_NUMBER_OF_ITERATIONS; ++iteration) {
            previousMetricValue = population[0]->partition->getMetricValue();
            applyCross(futures, parentIndices, population, denominator, step, amountOfNew);

            if ((iteration + 1) % config.MUTATION_FREQUENCY == 0) {
                ///get number of individuals for mutations
                unsigned long long amountOfMutants = Utils::generateRandomNumber(0,
                                                                                 (unsigned long long)(config.PERCENTAGE_OF_MUTANTS * amountOfNew));
                std::cout << "Amount of mutants: " << amountOfMutants << ' ';

                for (auto& mutation : mutations) {
                    ///shuffle for random selection of individuals for mutations
                    std::shuffle(population.begin() + sizeOfPopulation, population.end(), generator);

                    ///mutate
                    for (int individualIndex = 0; individualIndex < amountOfMutants; ++individualIndex) {
                        mutation(*population[individualIndex + sizeOfPopulation], config);
                    }
                }
            }

            if ((iteration + 1) % config.REFINE_FREQUENCY == 0) {
                ///greed optimization of partitions
                refinePopulation(futures, population, sizeOfPopulation, config);
            }

            ///sort by metric on partition to select a few best
            std::sort(population.begin(), population.end(),
                      [](const pGraph& left, const pGraph& right)
                      { return left->partition->getMetricValue() > right->partition->getMetricValue();});

            std::cout << "Score: " << std::setprecision(15) << population[0]->partition->getMetricValue() << '\n';

            ///algorithm is stopped when desired accuracy is got
            if (abs(previousMetricValue - population[0]->partition->getMetricValue()) < config.EPSILON)
                ++current_schedule;
            else
                current_schedule = 0;

            if (current_schedule == 0 && config.SAVE_INTERMEDIATE == 1) {
                std::string outputFileName(config.PREFIX + std::to_string(iteration) + ".out");
                population[0]->writePartition(outputFileName);
            }

            if (current_schedule == config.SCHEDULE) break;
        }

        graph.partition = population[0]->partition;
        graph.nodes = std::move(population[0]->nodes);
    }
}