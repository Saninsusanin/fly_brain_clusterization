#ifndef HGMC_COMMAND_LINE_H
#define HGMC_COMMAND_LINE_H

#include <string>
#include <random>

enum class Algorithm {
    LOUVAIN,
    LOUVAIN_STOCHASTIC,
    LEIDEN,
    K_NEAREST,
    GENETIC,
    RANDOM,
    NONE
};

enum class Mutation {
    MERGING,
    EXTRACTING,
    SEPARATING
};

enum class Engine {
    RANDOM_SPLIT,
    GIRVAN_NEWMAN
};

struct InputParameters {
    std::string inputFile;
    bool isMultiGraph;
    std::string outputFile;
    std::string configFile;
    Algorithm algorithm;

    Algorithm populationAlgorithm;
    std::vector<Mutation> mutation;
    unsigned int seed;
    std::string logFile;
    int verbose;
    Engine engine;

    explicit InputParameters(std::string inputFile=std::string(),
                             std::string outputFile=std::string(),
                             bool isMultiGraph=false,
                             std::string logFile=std::string(),
                             int verbose=0,
                             Algorithm algorithm=Algorithm::NONE,
                             Algorithm populationAlgorithm=Algorithm::K_NEAREST,
                             std::vector<Mutation> mutation=std::vector<Mutation> {},
                             unsigned int seed=std::random_device()(),
                             Engine engine=Engine::RANDOM_SPLIT);
};

InputParameters parseCommandLineArgs(int argc, char** argv);

#endif //HGMC_COMMAND_LINE_H
