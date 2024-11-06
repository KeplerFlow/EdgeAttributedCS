#ifndef ONLINE_H
#define ONLINE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <climits>
#include <queue>
#include <set>
#include <omp.h> 
#include "../Graph.h"

// 删除一个节点,检查一次连通性?
// 到达最小度数后,再删除所有不连通的节点
class CommunitySearch
{
public:
    static std::unordered_set<int> online(Graph &originalGraph, std::vector<int> &queryNodes, int k);

    static std::unordered_set<int> onlineWithEdge(Graph &originalGraph, std::vector<int> &queryNodes, int k,std::unordered_set<std::string> &lables);

    static std::unordered_set<int> onlineWithEdge(Graph &originalGraph, int k, std::unordered_set<std::string> &lables);

    static std::unordered_set<int> onlineWithCore(Graph &tempgraph, int k, std::unordered_set<std::string> &queryLabels,std::vector<std::unordered_set<int>> &core);

    static bool isConnected(Graph &graph, const std::vector<int> &queryNodes, const std::unordered_set<int> &subgraphNodes);
};

#endif