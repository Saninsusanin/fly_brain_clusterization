#include "utils.h"

#include <random>
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace kv;

std::mt19937 generator(12);

unsigned int Utils::generateRandomNumber(unsigned int min, unsigned int max) {
    std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);

    return dist(generator);
}


void Utils::setGenerator(unsigned int seed) {
   generator = std::mt19937(seed);
}


static int findInt(std::string& line) {
    auto pos = line.find_first_of("0123456789");
    if (pos == std::string::npos)
        throw std::invalid_argument("Incorrect genetic config");
    return std::stoi(line.substr(pos));
}


static double findDouble(std::string& line) {
    auto pos = std::min(line.find_first_of("123456789"), line.find_first_of('.'));
    if (pos == std::string::npos)
        throw std::invalid_argument("Incorrect genetic config");
    return std::stod(line.substr(pos));
}


GeneticConfig Utils::readGeneticConfig(std::string& pathToFile) {
    std::string line;
    std::fstream fin(pathToFile, std::ios::in);

    if (!fin.is_open())
        throw std::invalid_argument("No file with genetic config");

    int minSizeOfPopulation = -1;
    int maxSizeOfPopulation = -1;
    int maxNumberOfIterations = -1;
    int numberOfThreads = -1;
    double percentageOfMutants = -1.0;
    double epsilon = -1.0;
    int schedule = -1;
    int refineFrequency = 1;
    int mutationFrequency = 1;
    int minMergingCommunities = 2;
    double maxMergingPart = 0.25;
    int minExtractingCommunities = 1;
    double maxExtractingCommunitiesPart = 0.25;
    int minExtractingNodes = 1;
    double maxExtractingNodesPart = 0.25;
    int numberOfClusters = 10;
    int minSeparatingCommunities = 1;
    double maxSeparatingCommunitiesPart = 0.25;
    int saveIntermediate = 1;
    std::string prefix = "./";

    while (std::getline(fin, line)) {
        if (!line.empty()) {
            if (line[0] == '#')
                continue;
        }
        boost::algorithm::to_lower(line);

        if (line.find("min_size_of_population") != std::string::npos) {
            minSizeOfPopulation = findInt(line);
        }
        if (line.find("max_size_of_population") != std::string::npos) {
            maxSizeOfPopulation = findInt(line);
        }
        if (line.find("max_number_of_iterations") != std::string::npos) {
            maxNumberOfIterations = findInt(line);
        }
        if (line.find("number_of_threads") != std::string::npos) {
            numberOfThreads = findInt(line);
        }
        if (line.find("percentage_of_mutants") != std::string::npos) {
            percentageOfMutants = findDouble(line);
        }
        if (line.find("epsilon") != std::string::npos) {
            epsilon = findDouble(line);
        }
        if (line.find("schedule") != std::string::npos) {
            schedule = findInt(line);
        }
        if (line.find("refine_frequency") != std::string::npos) {
            refineFrequency = findInt(line);
        }
        if (line.find("mutation_frequency") != std::string::npos) {
            mutationFrequency = findInt(line);
        }
        if (line.find("min_merging_communities") != std::string::npos) {
            minMergingCommunities = findInt(line);
        }
        if (line.find("max_merging_part") != std::string::npos) {
            maxMergingPart = findDouble(line);
        }
        if (line.find("min_extracting_communities") != std::string::npos) {
            minExtractingCommunities = findInt(line);
        }
        if (line.find("max_extracting_communities_part") != std::string::npos) {
            maxExtractingCommunitiesPart = findDouble(line);
        }
        if (line.find("min_extracting_nodes") != std::string::npos) {
            minExtractingNodes = findInt(line);
        }
        if (line.find("max_extracting_nodes_part") != std::string::npos) {
            maxExtractingNodesPart = findDouble(line);
        }
        if (line.find("number_of_clusters") != std::string::npos) {
            numberOfClusters = findInt(line);
        }
        if (line.find("min_separating_communities") != std::string::npos) {
            minSeparatingCommunities = findInt(line);
        }
        if (line.find("max_separating_communities_part") != std::string::npos) {
            maxSeparatingCommunitiesPart = findDouble(line);
        }
        if (line.find("save_intermediate") != std::string::npos) {
            saveIntermediate = findInt(line);
        }
        if (line.find("prefix") != std::string::npos) {
            auto pos = line.find_first_of('=');
            prefix = line.substr(pos + 2);
        }
    }
    if (minSizeOfPopulation < 0 || maxSizeOfPopulation < 0 ||
            maxNumberOfIterations < 0 || numberOfThreads < 0 ||
            percentageOfMutants < 0. || percentageOfMutants > 1.0 ||
            epsilon < 0 || schedule < 0 || refineFrequency < 0 ||
            mutationFrequency < 0 || minMergingCommunities < 2 ||
            maxMergingPart < 0.0 || maxMergingPart > 1.0 ||
            minExtractingCommunities < 1 || maxExtractingCommunitiesPart < 0.0 ||
            maxExtractingCommunitiesPart > 1.0 || minExtractingNodes < 0 ||
            maxExtractingNodesPart < 0.0 || maxExtractingNodesPart > 0.5 ||
            numberOfClusters < 1 || minSeparatingCommunities < 1 ||
            maxSeparatingCommunitiesPart > 1.0)
        throw std::invalid_argument("Incorrect genetic config");

    return GeneticConfig(minSizeOfPopulation, maxSizeOfPopulation,
                         maxNumberOfIterations, numberOfThreads,
                         percentageOfMutants, epsilon, schedule,
                         refineFrequency, mutationFrequency,
                         minMergingCommunities, maxMergingPart,
                         minExtractingCommunities, maxExtractingCommunitiesPart,
                         minExtractingNodes, maxExtractingNodesPart,
                         numberOfClusters, minSeparatingCommunities,
                         maxSeparatingCommunitiesPart, saveIntermediate, prefix);
}

