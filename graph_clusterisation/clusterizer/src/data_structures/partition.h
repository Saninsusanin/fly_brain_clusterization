#ifndef HGMC_PARTITION_H
#define HGMC_PARTITION_H

#include "nodes.h"
#include "community.h"

namespace kv {
    struct Partition;
    using pPartition = std::shared_ptr<Partition>;

    struct Partition {
        double modularity;
        double regularization;
        std::vector<pCommunity> partition;

        explicit Partition(double modularity=0., double regularization=0.,
                           std::vector<pCommunity> partition=std::vector<pCommunity>());
        ~Partition() = default;

        void setMetricValue(int vertices, int edges) noexcept;
        double getMetricValue() const noexcept;
        void setModularity() noexcept;
        void updatePartition(Graph& graph) noexcept;
    };
}

#endif //HGMC_PARTITION_H
