#include <iostream>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "metric.h"


namespace po = boost::program_options;

struct Arguments {
    std::string realPartitionFile;
    std::string fakePartitionFile;
    std::string outputFile;

    explicit Arguments(std::string realPartitionFile = std::string(),
                       std::string fakePartitionFile = std::string(),
                       std::string outputFile = std::string()) {
        this->realPartitionFile = std::move(realPartitionFile);
        this->fakePartitionFile = std::move(fakePartitionFile);
        this->outputFile = std::move(outputFile);
    }
};


Arguments parseArguments(int argc, char **argv) {
    po::options_description desc("General options");
    std::string realPartitionFile;
    std::string fakePartitionFile;
    std::string outputFile("metric_value.out");

    desc.add_options()
            ("help,h", "Show help")
            ("real_partition,r", po::value<std::string>(&realPartitionFile), "Real labels")
            ("fake_partition,f", po::value<std::string>(&fakePartitionFile), "Predicted labels")
            ("output,o", po::value<std::string>(&outputFile),
             "Path to the file which will contain metric value");

    po::variables_map vm;

    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
        po::store(parsed, vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return Arguments();
        }
        if (realPartitionFile.empty()) {
            throw std::invalid_argument("No file with real partition was provided");
        }
        if (fakePartitionFile.empty()) {
            throw std::invalid_argument("No file with predicted partition was provided");
        }
    } catch (std::exception& exception) {
        std::cout << "Incorrect command" << std::endl;
    }
    return Arguments(realPartitionFile, fakePartitionFile, outputFile);
}


/**
 * example of real data path - "data/mel_2/real_l"
 * example of predicted data path - "../graph_clusterisation/clusterizer/output.out"
 * **/


int main(int argc, char* argv[]) {
    auto arguments = parseArguments(argc, argv);
    auto real_partition = read_partition(arguments.realPartitionFile);
    auto fake_partition = read_partition(arguments.fakePartitionFile);

    std::fstream output_file;

    output_file.open(arguments.outputFile, std::ios::out);

    output_file << std::setprecision(15) << kv_metric(real_partition, fake_partition);

    return 0;
}