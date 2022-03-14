#ifndef HGMC_UTILS_H
#define HGMC_UTILS_H

#include "../data_structures/nodes.h"
#include "../data_structures/genetic_config.h"
#include <random>

namespace kv {
    struct Utils {
        static unsigned int generateRandomNumber(unsigned int min, unsigned int max);
        static GeneticConfig readGeneticConfig(std::string& pathToFile);
        static void setGenerator(unsigned int seed);
    };
}

#endif //HGMC_UTILS_H