#include "graph.h"
#include <iostream>

namespace kv {
    Pair::Pair(int vertexId, int edges) {
        this->vertexId = vertexId;
        this->edges = edges;
    }


    Graph::Graph(std::vector<pNode> nodes, pAdjacencyList adjacencyList,
                 int edgesNumber, pPartition partition) {
        this->nodes = std::move(nodes);
        this->adjacencyList = std::move(adjacencyList);
        this->partition = std::move(partition);
        this->edgesNumber = edgesNumber;
    }

    Graph::Graph() {
        this->nodes = std::vector<pNode>();
        this->adjacencyList = nullptr;
        this->partition = nullptr;
        this->edgesNumber = 0;
    }


    int Graph::getVerticesNumber() const noexcept {
        return this->nodes.size();
    }

    const pNode& Graph::getNodeById(int nodeId) const noexcept {
        return this->nodes[nodeId - 1];
    }

    pGraph Graph::copy() const noexcept {
        std::vector<pNode> newNodes(this->nodes.size());

        for (int index = 0; index < newNodes.size(); index++) {
            newNodes[index] = this->nodes[index]->copy();
        }

        return std::make_shared<Graph>(std::move(newNodes), this->adjacencyList, this->edgesNumber);
    }

    const std::vector<Pair>& Graph::getNeighbours(int vertexId) const noexcept {
        return (*this->adjacencyList)[vertexId - 1];
    }

    void Graph::readGraphWrapped(std::istream &fin, Graph &graph, bool isMultiGraph) {
        if (!fin) {
            std::cout << "No such file has been found\n";
            return;
        }

        bool isFirstEdge = true;
        auto adjacencyList = std::make_shared<std::vector<std::vector<Pair>>>();
        int maxVertex = 0;

        while (!fin.eof()) {
            int vertex1, vertex2;
            int edges = 1;

            fin >> vertex1;
            fin >> vertex2;

            if (isMultiGraph)
                fin >> edges;

            if (isFirstEdge) {
                adjacencyList->resize((unsigned long)std::max(vertex1, vertex2) + 1);
                isFirstEdge = false;
            }

            if (std::max(vertex1, vertex2) >= adjacencyList->size()) {
                adjacencyList->resize((unsigned long)std::max(vertex1, vertex2) * 2);
            }

            (*adjacencyList)[vertex1].push_back(Pair(vertex2 + 1, edges));
            (*adjacencyList)[vertex2].push_back(Pair(vertex1 + 1, edges));
            maxVertex = std::max(maxVertex, std::max(vertex1, vertex2));
        }
        adjacencyList->resize(maxVertex + 1);

        std::vector<pNode> nodes(adjacencyList->size());
        int edgesNumber = 0;
        for (int index = 0; index < nodes.size(); index++) {
            nodes[index] = std::make_shared<Node>(index + 1,
                                                  0,
                                                  (int)(*adjacencyList)[index].size(),
                                                  0);
            edgesNumber += (*adjacencyList)[index].size();
        }

        graph.nodes = std::move(nodes);
        graph.adjacencyList = adjacencyList;
        graph.edgesNumber = edgesNumber / 2;
    }

    std::ostream &operator<<(std::ostream &fout, const Graph &graph) {
        if (graph.partition == nullptr) {
            return fout;
        }

        fout.precision(15);
        fout << graph.calculateMetricValue() << std::endl;
        for (const auto& community : graph.partition->partition) {
            for (const auto& node : community->community) {
                fout << node->vertexId - 1 << " ";
            }
            fout << std::endl;
        }

        return fout;
    }

    pGraph Graph::readGraph(std::string &inputFile, bool isMultiGraph) noexcept {
        auto graph = std::make_shared<Graph>();
        std::ifstream fin(inputFile);

        readGraphWrapped(fin, *graph, isMultiGraph);

        fin.close();
        return graph;
    }

    void Graph::writePartition(std::string &outputFile) const noexcept {
        std::ofstream fout(outputFile);

        fout << *this;

        fout.close();
    }

    double Graph::calculateDensity(const pCommunity &community) const noexcept {
        int size = (int)community->community.size();

        if (size == 1) {
            return 1.;
        }

        int innerEdges = community->getInnerEdges(*this);

        return 2. * innerEdges / size / (size - 1);
    }

    double Graph::calculateRegularization() const noexcept {
        double commonDensity = 0.;
        int vertices = this->getVerticesNumber();

        for (const auto& community : this->partition->partition) {
            commonDensity += calculateDensity(community);
        }

        int size = (int)partition->partition.size();

        return 0.5 * (commonDensity / size - double(size) / vertices);
    }

    double Graph::calculateMetricValue() const noexcept {
        double modularity = 0.;

        for (const auto& community : this->partition->partition) {
            int innerEdges = community->getInnerEdges(*this);
            int degree = community->getCommunityDegree();
            modularity += double(innerEdges) / this->edgesNumber -
                          pow(0.5 * degree / this->edgesNumber, 2);
        }
        //TODO:
        double reg = 0.;//calculateRegularization();
        return modularity + reg;
    }

    ///change nodes in this to louvain nodes
    void Graph::fromNodeGraphToLouvainGraph() noexcept {
        std::vector<pNode> newNodes(this->nodes.size());

        for (int index = 0; index < newNodes.size(); index++) {
            newNodes[index] = std::make_shared<LouvainNode>(this->nodes[index]);
        }

        this->nodes = std::move(newNodes);
    }
}