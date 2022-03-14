#ifndef HGMC_GRAPH_H
#define HGMC_GRAPH_H

#include "nodes.h"
#include "partition.h"
#include "community.h"

namespace kv {
    struct Graph;
    struct Pair;
    using pGraph = std::shared_ptr<Graph>;
    using pAdjacencyList = std::shared_ptr<std::vector<std::vector<Pair>>>;

    struct Pair{
        int vertexId;
        int edges;

        explicit Pair(int vertexId=0, int edges=0);
    };

    struct Graph {
        std::vector<pNode> nodes;
        pPartition partition;
        pAdjacencyList adjacencyList;
        int edgesNumber;

        Graph();
        Graph(std::vector<pNode> nodes, pAdjacencyList adjacencyList,
              int edgesNumber=0, pPartition partition=nullptr);

        int getVerticesNumber() const noexcept;
        const pNode& getNodeById(int nodeId) const noexcept;
        pGraph copy() const noexcept;
        const std::vector<Pair>& getNeighbours(int vertexId) const noexcept;

        friend std::ostream& operator<<(std::ostream& fout, const Graph& graph);

        static void readGraphWrapped(std::istream& fin, Graph& graph, bool isMultiGraph);
        static pGraph readGraph(std::string& inputFile, bool isMultiGraph=false) noexcept;
        void writePartition(std::string& outputFile) const noexcept;

        double calculateDensity(const pCommunity& community) const noexcept;
        double calculateRegularization() const noexcept;
        double calculateMetricValue() const noexcept;

        void fromNodeGraphToLouvainGraph() noexcept;
    };
}

#endif //HGMC_GRAPH_H
