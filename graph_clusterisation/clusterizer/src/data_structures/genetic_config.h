#ifndef HGMC_GENETIC_CONFIG_H
#define HGMC_GENETIC_CONFIG_H


struct GeneticConfig {
    int MIN_SIZE_OF_POPULATION;
    int MAX_SIZE_OF_POPULATION;
    int MAX_NUMBER_OF_ITERATIONS;
    int NUMBER_OF_THREADS;
    double PERCENTAGE_OF_MUTANTS;
    double EPSILON;
    int SCHEDULE;
    int REFINE_FREQUENCY;
    int MUTATION_FREQUENCY;
    int MIN_MERGING_COMMUNITIES;
    double MAX_MERGING_PART;
    int MIN_EXTRACTING_COMMUNITIES;
    double MAX_EXTRACTING_COMMUNITIES_PART;
    int MIN_EXTRACTING_NODES;
    double MAX_EXTRACTING_NODES_PART;
    int NUMBER_OF_CLUSTERS;
    int MIN_SEPARATING_COMMUNITIES;
    double MAX_SEPARATING_COMMUNITIES_PART;
    int SAVE_INTERMEDIATE;
    std::string PREFIX;
    
    explicit GeneticConfig(int minSizeOfPopulation, int maxSizeOfPopulation,
                           int maxNumberOfIterations, int numberOfThreads,
                           double percentageOfMutants, double epsilon, int schedule,
                           int refineFrequency, int mutationFrequency, int minMergingCommunities,
                           double maxMergingPart, int minExtractingCommunities, double maxExtractingCommunitiesPart,
                           int minExtractingNodes, double maxExtractingNodesPart, int numberOfClusters,
                           int minSeparatingCommunities, double maxSeparatingCommunitiesPart,
                           int saveIntermediate, std::string prefix) {
        MIN_SIZE_OF_POPULATION = minSizeOfPopulation;
        MAX_SIZE_OF_POPULATION = maxSizeOfPopulation;
        MAX_NUMBER_OF_ITERATIONS = maxNumberOfIterations;
        NUMBER_OF_THREADS = numberOfThreads;
        PERCENTAGE_OF_MUTANTS = percentageOfMutants;
        EPSILON = epsilon;
        SCHEDULE = schedule;
        REFINE_FREQUENCY = refineFrequency;
        MUTATION_FREQUENCY = mutationFrequency;
        MIN_MERGING_COMMUNITIES = minMergingCommunities;
        MAX_MERGING_PART = maxMergingPart;
        MIN_EXTRACTING_COMMUNITIES = minExtractingCommunities;
        MAX_EXTRACTING_COMMUNITIES_PART = maxExtractingCommunitiesPart;
        MIN_EXTRACTING_NODES = minExtractingNodes;
        MAX_EXTRACTING_NODES_PART = maxExtractingNodesPart;
        NUMBER_OF_CLUSTERS = numberOfClusters;
        MIN_SEPARATING_COMMUNITIES = minSeparatingCommunities;
        MAX_SEPARATING_COMMUNITIES_PART = maxSeparatingCommunitiesPart;
        SAVE_INTERMEDIATE = saveIntermediate;
        PREFIX = prefix;
    }
};


#endif