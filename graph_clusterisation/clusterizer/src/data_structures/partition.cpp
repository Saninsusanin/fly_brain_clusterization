#include "partition.h"
#include <utility>
#include <vector>

namespace kv {
    void Partition::setMetricValue(int vertices, int edges) noexcept {
        ///calculate modularity
        this->modularity = 0.;

        for (const auto& community : this->partition) {
            this->modularity += Community::calculateCommunityModularity(edges, community->innerEdges, community->communityDegree);
        }

        /// calculate regularization
        this->regularization = 0.;
    }

    double Partition::getMetricValue() const noexcept {
        return this->modularity + this->regularization;
    }

    void Partition::setModularity() noexcept {
        this->modularity = 0.;

        for (auto& community : this->partition) {
            this->modularity += community->modularity;
        }
    }


    void Partition::updatePartition(Graph& graph) noexcept {
        for (auto&& currentCommunity : this->partition) {
            currentCommunity->updateInnerEdges(graph);
            currentCommunity->updateCommunityDegree();
        }
    }

    Partition::Partition(double modularity, double regularization, std::vector<pCommunity> partition) {
        this->modularity = modularity;
        this->regularization = 0.;
        this->partition = std::move(partition);
    }

}