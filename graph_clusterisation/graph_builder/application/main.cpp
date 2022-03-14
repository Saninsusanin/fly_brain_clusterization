#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "graph_builder.h"
#include "distances.h"

namespace po = boost::program_options;

struct Arguments {
    std::string inputFile;
    std::string outputFile;
    std::string config;
    std::function<double(std::vector<double>&, std::vector<double>&)> distanceFunction;

    explicit Arguments(std::string inputFile = std::string(),
                       std::string outputFile = std::string(),
                       std::string configFile = std::string(),
                       std::function<double(std::vector<double>&, std::vector<double>&)> distanceFunction =
                               [](std::vector<double>& first_point, std::vector<double>& second_point) -> double {
        return get_minkovskiy_distance(first_point, second_point, 2); })  {
        this->inputFile = std::move(inputFile);
        this->outputFile = std::move(outputFile);
        this->config = std::move(configFile);
        this->distanceFunction = std::move(distanceFunction);
    };
};


Arguments parseArguments(int argc, char **argv) {
    po::options_description desc("General options");
    std::string inputFile;
    std::string outputFile("graph.out");
    std::string config;
    std::string distance;
    std::string p;
    double gamma;
    double bias;
    unsigned int power;
    std::function<double(std::vector<double>&, std::vector<double>&)> distanceFunction;

    desc.add_options()
            ("help,h", "Show help")
            ("input,i", po::value<std::string>(&inputFile), "Input file")
            ("output,o", po::value<std::string>(&outputFile), "Output file")
            ("config,c", po::value<std::string>(&config), "Config file: mean size of cluster,\n"
                                                          "weights number, weights")
            ("distance,dist", po::value<std::string>(&distance), "There exists several values for this parameter:\n"
                                                              "minkowski, polynomial, sigmoidal, rbf, laplacian")
            ("p", po::value<std::string>(&p), "Real number(excluding zero) or inf.\n"
                                              "Essential for minkovski distance")
            ("gamma", po::value<double>(&gamma), "Real number. Crucial for polynomial,\n"
                                                 "rbf, laplacian and sigmoidal kernels")
            ("bias", po::value<double>(&bias), "Real number. Crucial for polynomial and\n"
                                               "sigmoidal kernels")
            ("power", po::value<unsigned int>(&power), "Natural number.\n"
                                                       "Essential for polynomial kernel")
            ;

    po::variables_map vm;

    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
        po::store(parsed, vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return Arguments();
        } if (inputFile.empty()) {
            throw std::invalid_argument("No input file was provided");
        } if (config.empty()) {
            throw std::invalid_argument("No config file was provided");
        } if (distance.empty()) {
            throw std::invalid_argument("Empty distance parameter");
        } if (distance == "minkowski") {
            if (p == "inf") {
                distanceFunction = [&p](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return get_minkovskiy_distance(first_point, second_point, p);};
            } else {
                distanceFunction = [p](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return get_minkovskiy_distance(first_point, second_point, std::stod(p));};
            }
        } if (distance == "polynomial") {
            distanceFunction = [gamma, bias, power](std::vector<double>& first_point, std::vector<double>& second_point) {
                std::function<double(std::vector<double>&, std::vector<double>&)> kernel = [gamma, bias, power](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return polynomial_kernel(first_point, second_point, gamma, bias, power);
                };
                return get_kernel_distance(first_point, second_point, kernel);
            };
        } if (distance == "sigmoidal") {
            distanceFunction = [gamma, bias](std::vector<double>& first_point, std::vector<double>& second_point) {
                std::function<double(std::vector<double>&, std::vector<double>&)> kernel = [gamma, bias](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return sigmoid_kernel(first_point, second_point, gamma, bias);
                };
                return get_kernel_distance(first_point, second_point, kernel);
            };
        } if (distance == "rbf") {
            distanceFunction = [gamma](std::vector<double>& first_point, std::vector<double>& second_point) {
                std::function<double(std::vector<double>&, std::vector<double>&)> kernel = [gamma](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return exponential_kernel(first_point, second_point, gamma, 2);
                };
                return get_kernel_distance(first_point, second_point, kernel);
            };
        } if (distance == "laplacian") {
            distanceFunction = [gamma](std::vector<double>& first_point, std::vector<double>& second_point) {
                std::function<double(std::vector<double>&, std::vector<double>&)> kernel = [gamma](std::vector<double>& first_point, std::vector<double>& second_point) {
                    return exponential_kernel(first_point, second_point, gamma, 1);
                };
                return get_kernel_distance(first_point, second_point, kernel);
            };
        }
    } catch (std::exception& exception) {
        std::cout << "Incorrect command" << std::endl;
    }
    return Arguments(inputFile, outputFile, config, distanceFunction);
}

struct Params {
    std::vector<int> weights;
    double size;

    explicit Params(std::vector<int> weights, double size) {
        this->weights = std::move(weights);
        this->size = size;
    }
};

Params readConfig(std::string& config) {
    double size;
    int weightsNum;
    std::vector<int> weights;
    
    std::ifstream fin(config);
    fin >> size;
    fin >> weightsNum;
    for (int i = 0; i < weightsNum; i++) {
        int tmp = 0;
        fin >> tmp;
        weights.push_back(tmp);
    }
    return Params(weights, size);
}


int main(int argc, char** argv) {
    auto arguments = parseArguments(argc, argv);
    auto params = readConfig(arguments.config);

    auto graph = get_graph(arguments.inputFile, params.size, params.weights, arguments.distanceFunction);
    print_graph(graph, arguments.outputFile);

    return 0;
}
