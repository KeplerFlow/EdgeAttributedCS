#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include <map>
#include <queue>
class Graph
{
public:
    // 构造函数和析构函数

    // 从文件中读取图数据，并计算节点的度和最小度
    Graph(std::string path);
    Graph(std::string path,int radio);
    Graph(std::string path,int v,int y,std::string label);
    // 构造一个拷贝图
    Graph(Graph &graph);
    ~Graph(); // 析构函数

    // 基本图操作
    bool addEdge(int from, int to); // 添加两个节点之间的边

    void addNode(int node); // 添加节点到图中

    std::unordered_set<std::string> getEdgeLabels(int v,int n);

    std::unordered_map<int, int> computeDegrees(); // 计算每个节点的度

    int computeMinimumDegree(); // 计算图中的最小度

    std::vector<std::pair<int, int>> getEdges(); // 返回图的边

    bool doesNodeExist(int node); // 检查节点是否存在于图中

    std::unordered_map<int, int> getDegrees(); // 获取所有节点的度

    int getDegree(int node); // 获取特定节点的度

    int getMinimumDegree(); // 获取图中的最小度

    std::unordered_set<int> getNeighbors(int node); // 获取特定节点的邻居节点

    std::vector<int> sortNeighbors(int node);

    int getNumberOfNodes(); // 获取图中节点的数量

    int getNumberOfEdges(); // 获取图中边的数量

    int getNumberOfEdgesOfConnectedComponent(std::unordered_set<int> connectedComponent); // 获取连接组件中边的数量

    std::vector<std::unordered_set<int>> getOrderedNodes();

    std::unordered_set<int> getNodes(); // 获取图中所有节点

    void removeNode(int node);

    void removeEdge(int from, int to, std::string lable);

    void serializeToFile(std::string datasetName); // 将图序列化到文件中

    void printDegrees();

    void printOrderedNodes();

    void filteEdges(std::string label);

    std::unordered_set<int> findLargestConnectedComponent();

    void printGraphInfo(Graph &graph);

    bool isGraphConnected();

    std::vector<std::unordered_set<int>> computeConnectedComponents();

    void keepLargestComponent();

    std::unordered_set<std::string> getNodesLable(int id);

    std::unordered_set<std::string> getEdgeLable(int from, int to);

    void addEdgeLable(int from, int to, std::string lable);

    void addNodeLable(int id, std::string lable);

    void filterEdgeLable(std::unordered_set<std::string> lables);

    void filterEdgeLableAndRecover(std::unordered_set<std::string> lables);

    std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>>& getEdgeLables();
    
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> getLabelsForEdges();

    std::vector<std::pair<int, int>> addEdgeSet(std::vector<std::pair<int, int>> edgeSet,std::string label);

    std::unordered_map<int, int> computingOrder();

    void statistic();

    int getNodeNum();

    std::unordered_set<std::string> checkLabelForm();

    int getEdgeNum();

    std::unordered_map<std::string, std::vector<std::pair<int, int>>>& getLabelsForDirectedEdges();

    bool isExistEdge(int v1,int v2);

    int getDegree();

    std::unordered_map<int, std::unordered_set<int>>& getAdj();

private:
    void dfsStack(int startNode, std::unordered_set<int> &visited);

    void dfs(int node, std::unordered_set<int> &visited, std::unordered_set<int> &component);

private:
    std::unordered_map<int, std::unordered_set<int>> adj; // 图的邻接表表示：节点 ID 和其相邻节点

    std::unordered_map<int, int> degrees; // 图中节点的度

    int Dmax;

    int m;

    int n;

    // 一个节点对应多个属性
    std::unordered_map<int, std::unordered_set<std::string>> nodeLables; // 图中节点的名称

    // 两个节点之间可能存在多边
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> edgeLables;

    std::vector<std::unordered_set<int>> orderedNodes; // 表示具有相同度的节点的集合的向量,key为degree

    int minimumDegree; // 图中所有节点的最小度

    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForEdges;//最外层键为标签,内层为该图中边

    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForDirectedEdges;//最外层键为标签,内层为该图中边,小的在前

    std::vector<std::vector<int>> O;//顺序存储每个Ok的节点

    std::unordered_map<int,std::pair<int,int>> OIndex;//每个节点在O中的索引

    // 辅助函数
    void readFromFile(std::string fileName); // 从文件中读取图
    void readFromFile(std::string fileName,int radio); // 从文件中读取图
    void readFromFile(std::string fileName,int v,int y,std::string label); // 从文件中读取图
};

#endif // GRAPH_H
