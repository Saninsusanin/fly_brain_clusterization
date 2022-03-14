#ifndef HGMC_COMMUNITY_H
#define HGMC_COMMUNITY_H

#include "nodes.h"

namespace kv {
    struct Graph;
    struct Community;
    struct CommunityComparator;

    using pCommunity = std::shared_ptr<Community>;
    using CommunityStorage = std::set<pNode, CommunityComparator>;


    struct CommunityComparator {
        bool operator() (const pNode& left, const pNode& right) const {
            return left->vertexId < right->vertexId;
        }
    };


    struct Community {
        double modularity;
        int communityDegree;
        int innerEdges;
        int aggregationNumber;
        CommunityStorage community;

        explicit Community(CommunityStorage community=CommunityStorage(),
                           int communityDegree=0, int innerEdges=0, double modularity=0., int aggregatedNumber=0);

        int getCommunityId() const noexcept;
        void updateModularity(int edges) noexcept;
        void updateInnerEdges(Graph& graph) noexcept;
        void updateCommunityDegree() noexcept;
        double getCommunityDensity() const noexcept;
        void updateAggregationNumber() noexcept;

        void insertNode(const pNode& node,
                        int innerEdgesFromNode,
                        int edges, int nodeAggregationNumber=1);
        void eraseNode(const pNode& node,
                       int innerEdgesFromNode,
                       int edges, int nodeAggregationNumber=1);

        int getInnerEdges(const Graph& graph) const noexcept;
        int getCommunityDegree() const noexcept;
        static double calculateCommunityModularity(int edgesNumber,
                                                   int _innerEdges,
                                                   int _communityDegree) noexcept;
    };
}

#endif //HGMC_COMMUNITY_H
