#include "algorithms/genetic_algorithm.h"
#include "algorithms/louvain_algorithm.h"
#include "algorithms/leiden_algorithm.h"
#include "data_structures/graph.h"
#include "command_line.h"
#include "utils/utils.h"
#include <chrono>
#include <functional>
#include <random>
#include <iostream>
#include "logger.h"

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

extern std::mt19937 generator;


void runLouvainLeiden(kv::pGraph& graph,
                      std::string& outputFile,
                      Algorithm& algorithm) {
    graph->fromNodeGraphToLouvainGraph();

    if (algorithm == Algorithm::LOUVAIN) {
        Louvain(*graph, true, true);
    }
    else if (algorithm == Algorithm::LOUVAIN_STOCHASTIC) {
        Louvain(*graph, false, true);
    }
    else if (algorithm == Algorithm::LEIDEN) {
        Leiden(*graph);
    }

    graph->writePartition(outputFile);
}

void runKNearest(kv::pGraph& graph, std::string& outputFile) {
    kNearestNeighboursPartition(*graph);

    graph->writePartition(outputFile);
}

void runGenetic(kv::pGraph& graph,
                std::string& outputFile,
                Algorithm& populationAlgorithm,
                std::vector<Mutation>& mutations,
                std::string& configFile,
                Engine& engine) {

    std::vector<std::function<void(kv::Graph&, GeneticConfig&)>> mutationMethods;

    auto geneticConfig = kv::Utils::readGeneticConfig(configFile);

    std::function<std::vector<std::vector<int>>(std::vector<kv::pNode>&)> engineFunction;

    if (engine == Engine::GIRVAN_NEWMAN) {
        auto& tmp = *graph;
        engineFunction = [&tmp](std::vector<kv::pNode>& community) {return kv::GirvanNewman(tmp, community);};
    } else if (engine == Engine::RANDOM_SPLIT) {
        engineFunction = kv::randomSplit;
    }

    for (auto mutation : mutations) {
        if (mutation == Mutation::MERGING) {
            mutationMethods.emplace_back(kv::mergingMutate);
        } else if (mutation == Mutation::EXTRACTING) {
            mutationMethods.emplace_back(kv::extractingMutate);
        } else if (mutation == Mutation::SEPARATING) {
            mutationMethods.emplace_back([&engineFunction](kv::Graph& graph, GeneticConfig& config) {
                return kv::separatingMutate(graph, config, engineFunction);
            });
        }
    }

    std::function<void(kv::Graph&)> populationGenerator;

    if (populationAlgorithm == Algorithm::LOUVAIN_STOCHASTIC) {
        populationGenerator = [](kv::Graph& graph) {return kv::Louvain(graph, false);};
    } else if (populationAlgorithm == Algorithm::K_NEAREST) {
        populationGenerator = kv::kNearestNeighboursPartition;
    } else if (populationAlgorithm == Algorithm::LEIDEN) {
        populationGenerator = kv::Leiden;
    } else if (populationAlgorithm == Algorithm::RANDOM) {
        populationGenerator = [&geneticConfig](kv::Graph& graph) {return kv::randomInitialization(graph, geneticConfig.NUMBER_OF_CLUSTERS);};
    }

    graph->fromNodeGraphToLouvainGraph();

    geneticAlgorithm(*graph, populationGenerator, mutationMethods, geneticConfig);

    graph->writePartition(outputFile);
}

void run(InputParameters& params) {
    if (params.algorithm == Algorithm::NONE) {
        return;
    }

    auto graph = kv::Graph::readGraph(params.inputFile, params.isMultiGraph);

    if (graph->adjacencyList == nullptr) {
        return;
    }

    if (params.algorithm == Algorithm::LOUVAIN ||
        params.algorithm == Algorithm::LEIDEN ||
        params.algorithm == Algorithm::LOUVAIN_STOCHASTIC) {
        runLouvainLeiden(graph, params.outputFile, params.algorithm);
    } else if (params.algorithm == Algorithm::K_NEAREST) {
        runKNearest(graph, params.outputFile);
    } else if (params.algorithm == Algorithm::GENETIC) {
        runGenetic(graph,
                   params.outputFile,
                   params.populationAlgorithm,
                   params.mutation,
                   params.configFile,
                   params.engine);
    }

}

int main(int argc, char** argv) {
    auto params = parseCommandLineArgs(argc, argv);
    // init(params.logFile, params.verbose);
    kv::Utils::setGenerator(params.seed);
    auto t1 = std::chrono::high_resolution_clock::now();
    run(params);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time in milliseconds: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << '\n';

    return 0;
}