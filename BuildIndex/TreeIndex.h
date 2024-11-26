#ifndef TREEINDEX_H
#define TREEINDEX_H

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

#include "../Graph.h"
#include "../Algorithm/CoreGroup.h"
class TreeIndex 
{
public:
    // 核心指标，记录每个节点所属的核心编号
    std::unordered_map<int, int> coreIndex;

    // 每个核心的最小度数
    std::unordered_map<int, int> coreMinimumDegree;

    // 第一个键为壳的层数,内层的键是连接组件ID，值是父连接组件的ID
    std::unordered_map<int, std::unordered_map<int, int>> connectedComponentParent;

    // 第一个键为壳的层数,第二个键为连接组件ID,set内存储子连接组件的ID
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> connectedComponentChildren;

    // 第一个键是壳的层数，内层的键是连接组件ID，值是一个HashSet，包含连接组件中的节点ID
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> connectedComponentNodes;

    // 节点分组，记录每个节点所属的connectedComponentNodes分组编号
    std::unordered_map<int, int> nodeGroup;

    // 键是节点的ID，内层的键是壳的核心指数，值是一个HashSet，包含节点在该壳中的邻居节点ID。
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> nodeNeighbors;

    std::unordered_map<int, std::unordered_map<int, std::unordered_set<int>>> layerToComponentToNodes;

    std::unordered_map<int, std::unordered_set<int>> ComponentToNodes;
    
    std::unordered_map<int, int> nodeToComponentId;

    std::unordered_map<int, std::unordered_set<int>> ComponentParent;

    std::unordered_map<int, std::unordered_set<int>> ComponentChildren;

    //<节点 <邻居父节点,标签 >>
    std::unordered_map<int,std::unordered_map<int,std::unordered_set<std::string>>> parentIndex;
    
    //<节点 <邻居子节点,标签 >>
    std::unordered_map<int,std::unordered_map<int,std::unordered_set<std::string>>> childIndex;
    
    int nextComponentId = 0; // 用于生成唯一的连通分量 ID

public:
    TreeIndex(Graph &graph,std::unordered_set<std::string> &labels);
    void buildParentChildRelationships(Graph &graph);
    void identifyAndStoreComponents(Graph &graph);
    void printComponents() const;
    void calculateParentRelationships(Graph &graph);
    void calculateChildRelationships(Graph &graph);
    std::unordered_set<int> findTopLevelComponents(std::vector<int> &queryNodes,int k);
    std::unordered_set<int> findTopComponents(std::vector<int> &queryNodes,int k);
    std::unordered_set<int> findMaxKCoreSubgraph(std::vector<int>& queryNodes, int k);
    std::unordered_set<int> findKCoreSubgraph(std::vector<int>& queryNodes, int k);
    void exploreDownwards(int compId, int currentLevel, int k, std::unordered_set<int>& resultNodes, std::queue<int>& componentQueue);

    // 构造函数，根据图构建树索引
    TreeIndex(Graph &graph, std::string datasetName);

    // 构造函数，从文件中加载树索引
    TreeIndex(const std::string &path);

    // 计算核心指标
    void computeCoreIndex(Graph &graph,std::unordered_set<std::string> &labels);

    // 计算核心组成
    void computeCoreComposition(Graph &graph);

    void computeCoreCompositionByLcy(Graph &graph);

    // 计算节点邻居信息
    void computeNodeNeighbors(Graph &graph);

    // 序列化树索引到文件
    void serializeToFile(const std::string &path);

    // 从文件中反序列化树索引
    void deserializeFromFile(const std::string &path);

    // 获取一组节点中的最小核心编号
    int getMinimumCoreIndex(std::vector<int> queryNodes);

    // 获取节点的邻居节点集合（在指定核心下）
    std::unordered_set<int> getNeighbors(int node, int coreIndex);

    // 获取节点的邻居节点集合（在指定核心下，并且在指定子核心中）
    std::unordered_set<int> getNeighbors(int node, int coreIndex, std::unordered_set<int> &subcore);

    // 获取指定核心的最小度数
    int getCoreMinimumDegree(int coreIndex);

    // 获取指定核心中的节点数量
    int getNumberOfNodes(int coreIndex);

    // 获取一组节点的核心集合
    std::unordered_set<int> getCore(std::vector<int> queryNodes);

    // 获取节点的核心编号
    int getCoreIndex(int node);
    // Get the component ID for a given node
    int getComponent(int node);

    // Get the parent component ID for a given component
    int getParentComponent(int componentId);
    std::unordered_set<int> getConnectedComponentChildren(int componentId);

    // Get all nodes in a specified component
    std::unordered_set<int> getNodesInComponent(int componentId);
    // 打印核心指标
    void printCoreIndex();

    // 打印每个核心的最小度数
    void printCoreMinimumDegree();

    // 打印连通分量的父子关系
    void printConnectedComponentParent();

    // 打印连通分量的孩子节点集合
    void printConnectedComponentChildren();

    // 打印连通分量的节点集合
    void printConnectedComponentNodes();
    // 打印节点分组
    void printNodeGroup();

    // 打印节点邻居信息
    void printNodeNeighbors();
    void printQueryNodesCoreIndex(const std::vector<int>& queryNodes);

};

#endif // TREE_H
