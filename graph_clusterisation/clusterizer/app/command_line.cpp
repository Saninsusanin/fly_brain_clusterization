#include "command_line.h"
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>
#include <string>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

namespace po = boost::program_options;

InputParameters::InputParameters(std::string inputFile,
                                 std::string outputFile,
                                 bool isMultiGraph,
                                 std::string logFile,
                                 int verbose,
                                 Algorithm algorithm,
                                 Algorithm populationAlgorithm,
                                 std::vector<Mutation> mutation,
                                 unsigned int seed,
                                 Engine engine) {
    this->inputFile = std::move(inputFile);
    this->outputFile = std::move(outputFile);
    this->isMultiGraph = isMultiGraph;
    this->algorithm = algorithm;
    this->populationAlgorithm = populationAlgorithm;
    this->mutation = mutation;
    this->configFile = std::string();
    this->seed = seed;
    this->logFile = std::move(logFile);
    this->verbose = verbose;
    this->engine = engine;
}

static InputParameters getLouvainParameters(std::string& inputFile,
                                            std::string& outputFile,
                                            bool isMultiGraph,
                                            std::string& logFile,
                                            int verbose,
                                            po::variables_map& vm) {
    InputParameters params(inputFile, outputFile, isMultiGraph, logFile, verbose);
    if (vm.count("stochastic")) {
        params.algorithm = Algorithm::LOUVAIN_STOCHASTIC;
    } else {
        params.algorithm = Algorithm::LOUVAIN;
    }

    if (vm.count("seed")) {
        params.seed = vm["seed"].as<unsigned int>();
    }

    return params;
}


static InputParameters getKNearestNeighboursParameters(std::string& inputFile,
                                                       std::string& outputFile,
                                                       bool isMultiGraph,
                                                       std::string& logFile,
                                                       int verbose,
                                                       po::variables_map& vm) {
    InputParameters params(inputFile, outputFile, isMultiGraph, logFile, verbose);

    if (vm.count("seed")) {
        params.seed = vm["seed"].as<unsigned int>();
    }

    return params;
}


static InputParameters getGeneticParameters(std::string& inputFile,
                                            std::string& outputFile,
                                            bool isMultiGraph,
                                            std::string& logFile,
                                            int verbose,
                                            po::variables_map& vm) {
    InputParameters params(inputFile, outputFile, isMultiGraph, logFile, verbose, Algorithm::GENETIC);

    if (!vm.count("config"))
        throw std::invalid_argument("No path to config");

    params.configFile = vm["config"].as<std::string>();

    if (vm.count("seed")) {
        params.seed = vm["seed"].as<unsigned int>();
    }

    if (vm.count("mut")) {
        auto mutationType = vm["mut"].as<std::vector<std::string>>();

        for (auto& mutation : mutationType) {
            if (mutation == "merging") {
                params.mutation.push_back(Mutation::MERGING);
            } else if (mutation == "extracting") {
                params.mutation.push_back(Mutation::EXTRACTING);
            } else if (mutation == "separating") {
                params.mutation.push_back(Mutation::SEPARATING);
            } else {
                params.algorithm = Algorithm::NONE;
                std::cout << "Incorrect mutation parameter" << std::endl;

                return params;
            }
        }
    }

    if (vm.count("engine")) {
        auto engine = vm["engine"].as<std::string>();

        if (engine == "random-split") {
            params.engine = Engine::RANDOM_SPLIT;
        } else if (engine == "girvan-newman") {
            params.engine = Engine::GIRVAN_NEWMAN;
        } else {
            std::cout << "Incorrect population generator algorithm" << std::endl;
            return params;
        }
    }

    if (vm.count("gen")) {
        auto populationAlgorithm = vm["gen"].as<std::string>();

        if (populationAlgorithm == "louvain") {
            params.populationAlgorithm = Algorithm::LOUVAIN_STOCHASTIC;
        } else if (populationAlgorithm == "leiden") {
            params.populationAlgorithm = Algorithm::LEIDEN;
        } else if (populationAlgorithm == "k-nearest") {
            params.populationAlgorithm = Algorithm::K_NEAREST;
        } else if (populationAlgorithm == "random") {
            params.populationAlgorithm = Algorithm::RANDOM;
        } else {
            params.algorithm = Algorithm::NONE;
            std::cout << "Incorrect population generator algorithm" << std::endl;
            return params;
        }
    }
    return params;
}

InputParameters parseCommandLineArgs(int argc, char** argv) {
    po::options_description desc("General options");
    std::string algo;
    std::string inputFile;
    std::string outputFile("output.out");
    std::string logFile("log.log");
    int verbose = 0;

    desc.add_options()
            ("help,h", "Show help")
            ("input,i", po::value<std::string>(&inputFile), "Input file")
            ("multi,m", "Does input file contain multigraph")
            ("output,o", po::value<std::string>(&outputFile), "Output file")
            ("algo,a", po::value<std::string>(&algo),
             "Select algorithm: louvain, leiden, k-nearest, genetic")
            ("log,l", po::value<std::string>(&logFile), "File with log")
            ("verbose,v", po::value<int>(&verbose), "logging level\n"
                                                    "6 - provide only score of the best partition of an each population of the genetic algorithm\n"
                                                    "5 - leave only parameters of the genetic algorithm\n"
                                                    "4 - leave only time measurement");

    po::options_description knn_desc("KNearestNeighbours options");
    knn_desc.add_options()
            ("seed", po::value<unsigned int>(), "seed for the random number engine");

    po::options_description louvain_desc("Louvain options");
    louvain_desc.add_options()
            ("stochastic", "Call stochastic variant")
            ("seed", po::value<unsigned int>(), "seed for the stochastic variant");

    std::string geneticConfigFile;
    po::options_description genetic_desc("Genetic options");
    genetic_desc.add_options()
            ("gen,g", po::value<std::string>(), "Population generator algorithm: "
                                                "louvain, leiden or k-nearest")
            ("seed", po::value<unsigned int>(), "seed for the random number engine")
            ("mut", po::value<std::vector<std::string>>()->multitoken(), "Mutation type: merging, extracting or separating\n"
                                                                         "you can pass several and they would be applied "
                                                                         "in the order in which they were mentiooned")
            ("config", po::value<std::string>(&geneticConfigFile), "Path to genetic config file")
            ("engine", po::value<std::string>(), "Backend for the separating mutate algorithm");

    po::variables_map vm;

    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
        po::store(parsed, vm);
        po::notify(vm);
        auto isMultiGraph = false;

        if (vm.count("help")) {
            std::cout << desc << '\n';
            std::cout << knn_desc << '\n';
            std::cout << louvain_desc << '\n';
            std::cout << genetic_desc << std::endl;

            return InputParameters(std::string(), std::string(), false, std::string(), 0, Algorithm::NONE);
        }
        if (inputFile.empty()) {
            std::cout << "No input file was provided" << std::endl;

            return InputParameters(std::string(), std::string(), false, std::string(), 0, Algorithm::NONE);
        }
        if (vm.count("multi")) {
            isMultiGraph = true;
        }
        if (vm.count("log")) {
            logFile = vm["log"].as<std::string>();
        }
        if (vm.count("verbose")) {
            verbose = 6 - vm["verbose"].as<int>();
        }
        if (algo == "louvain") {
            desc.add(louvain_desc);
            po::store(po::parse_command_line(argc, argv, desc), vm);

            return getLouvainParameters(inputFile, outputFile, isMultiGraph, logFile, verbose, vm);
        }
        else if (algo == "leiden") {
            return InputParameters(inputFile, outputFile, isMultiGraph, logFile, verbose, Algorithm::LEIDEN);
        }
        else if (algo == "k-nearest") {
            desc.add(knn_desc);
            po::store(po::parse_command_line(argc, argv, desc), vm);

            return getKNearestNeighboursParameters(inputFile, outputFile, isMultiGraph, logFile, verbose, vm);
        }
        else if (algo == "genetic") {
            desc.add(genetic_desc);
            po::store(po::parse_command_line(argc, argv, desc), vm);

            return getGeneticParameters(inputFile, outputFile, isMultiGraph,logFile, verbose, vm);
        } else {
            std::cout << "Incorrect algorithm type" << std::endl;

            return InputParameters(std::string(), std::string(), isMultiGraph, logFile, verbose, Algorithm::NONE);
        }
    } catch (std::exception& exception) {
        std::cout << "Incorrect command" << std::endl;
    }

    return InputParameters();
}