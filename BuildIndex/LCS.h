#ifndef LCS_H
#define LCS_H
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <climits>
#include <chrono>
#include <cassert>
#include <cstring>
#include <memory>
#include <cassert>

#include "../Graph.h"
#include "../Algorithm/CoreGroup.h"
#include "../core-maint/core.h"
#include "../core-maint/defs.h"
#include "../core-maint/gadget/gadget.h"
#include "../core-maint/glist/glist.h"
class LCS
{

private:
    // 最外层key为顶点id,次外层为k值,内层集合为边标签集,集合中存在多个不相交的边标签子集
    std::unordered_map<int, std::unordered_map<int, std::vector<std::unordered_set<std::string>>>> LcsIndex;
    // 边标签集枚举得到的非空子集,按照从少到多.
    std::vector<std::unordered_set<std::string>> subLabels;
    // 边标签集枚举得到的非空子集,按照标签集内元素个数.
    std::unordered_map<int, std::vector<std::unordered_set<std::string>>> bfsLabels;
    // key为标签集,即subLabels的元素,值为coreness键值对,
    std::unordered_map<int, std::unordered_map<int, int>> bfsCoreness;
    // 无序的边标签子集
    std::vector<std::unordered_set<std::string>> normalSubLabels;

    std::unordered_map<int,std::unordered_map<int, std::unordered_set<int>>> mupt_index;
    std::unordered_map<int,std::multimap<int, int, std::greater<>>> mupt_index_map;
    
    std::vector<int> indexs_dfs;
    std::unordered_map<int,std::unordered_set<std::string>> indexs_to_subset_dfs;
    std::unordered_map<int, int> treeIndexs_dfs;
    std::unordered_map<int, std::string> treeLabelIndexs_dfs;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForEdges_dfs;

    std::unordered_map<int, std::unordered_set<int>> subcore;
    std::unordered_map<int,int> cd;

private:
    void printSubLabels();
    bool isSubSet(std::vector<std::unordered_set<std::string>> &lcs, std::unordered_set<std::string> &labelSet);
    bool isSubsetOf(const std::unordered_set<std::string> &subset, const std::unordered_set<std::string> &set);

public:
    LCS();

    LCS(std::unordered_set<std::string>, Graph &graph);

    LCS(std::unordered_set<std::string>, Graph &graph, bool flag);

    std::unordered_map<int, std::unordered_map<int, std::vector<std::unordered_set<std::string>>>> getLcsIndex() { return LcsIndex; }

    void build(Graph &graph);

    void build_in_mupt_index(Graph &graph);

    void opt(Graph &graph);

    void optWithSpend(Graph &graph);

    void opt_dfs(Graph &graph);

    void dfs(Graph &graph,std::vector<std::vector<int>> graphForCore,std::vector<int> core,core::GLIST cm,std::unordered_set<std::string> &subset,int subIndex);

    // 枚举边标签集得到子集
    void enumLabels(std::unordered_set<std::string> labels);

    // 根据边标签集截断,得到每个顶点对应的coreness
    std::unordered_map<int, int> coreDecomposition(std::unordered_set<std::string> labels, Graph &graph);

    bool isExist(int vertex, std::unordered_set<std::string> &labels, int k);
    void printLcsIndex();

    std::unordered_set<int> coreMaintenance(Graph &graph,std::unordered_map<int, int> coreness,std::unordered_map<std::string, std::vector<std::pair<int, int>>> &labelsForEdges);

    void incIndex(Graph &graph,int v,int u,std::string label);
    void incIndex(std::unordered_map<int, std::unordered_set<int>> &adj,int v,int u,std::unordered_set<std::string> labels,std::string label);
    void findSubCore(std::unordered_map<int, std::unordered_set<int>> &adj,int u,std::unordered_set<std::string> labels);  
    int vertexToK(int u,std::unordered_set<std::string> labels);

    void serialize(std::ofstream &ofs) const ;
    void serialize_mupt_index(std::ofstream &ofs) const ;
    void deserialize(std::ifstream &ifs);
    void deserialize_mupt_index(std::ifstream &ifs);
    size_t calculateMemoryUsage();
    std::unordered_map<int,std::unordered_set<int>> mupt_index_search(std::unordered_set<std::string> labels);
    std::unordered_set<int> mupt_index_search_map(std::unordered_set<std::string> labels,int k);

    void transIndex();

};

#endif