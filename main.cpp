#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <string.h>
#include <unordered_set>
#include <random>
#include <ctime>
#include <fstream>
#include <filesystem> // C++17 for file existence check
#include "Graph.h"
#include "./Algorithm/CoreGroup.h"
#include "./BuildIndex/TreeIndex.h"
#include "./BuildIndex/hierarchy.h"
#include "./BuildIndex/LCS.h"
#include "./Algorithm/Online.h"
#include "./utils/Mapping.h"
#include "./core-maint/core.h"
#include "./core-maint/defs.h"
#include "./core-maint/gadget/gadget.h"
#include "./core-maint/glist/glist.h"
// 序列化索引到文件
void serializeIndex(const LCS &optLcs, const std::string &indexFile)
{
    std::ofstream ofs(indexFile, std::ios::binary);
    if (!ofs)
    {
        std::cerr << "Error opening file for writing: " << indexFile << std::endl;
        return;
    }
    // 这里假设LCS类可以被序列化（可以自己实现LCS的序列化函数）
    optLcs.serialize(ofs);
    ofs.close();
}

// 从文件反序列化读取索引
bool deserializeIndex(LCS &optLcs, const std::string &indexFile)
{
    std::ifstream ifs(indexFile, std::ios::binary);
    if (!ifs)
    {
        std::cerr << "Error opening file for reading: " << indexFile << std::endl;
        return false;
    }
    // 这里假设LCS类有一个反序列化的函数
    optLcs.deserialize(ifs);
    ifs.close();
    return true;
}

std::unordered_set<int> computeCore(Graph &graph, int querynode, int k)
{
    std::cout << "arrive";
    int n = graph.getNodeNum() + 1;
    std::vector<std::vector<int>> graphForCore;
    graphForCore.resize(n);
    std::vector<int> core;
    core.resize(n);
    auto cm = std::make_shared<core::GLIST>(n);
    std::vector<std::pair<int, int>> edges = graph.getEdges();

    for (int j = 0; j < edges.size(); ++j)
    {
        int v1 = edges.at(j).first;
        int v2 = edges.at(j).second;
        if (0 <= v1 && v1 < n && 0 <= v2 && v2 < n && v1 != v2)
        {
            graphForCore[v1].push_back(v2);
            graphForCore[v2].push_back(v1);
        }
    }
    cm->ComputeCore(graphForCore, true, core);

    if (core[querynode] < k)
    {
        return {};
    }

    // 从查询节点开始检索最大连通分量
    std::queue<int> q;
    std::unordered_set<int> visited;
    std::unordered_set<int> largestComponent;

    q.push(querynode);
    visited.insert(querynode);

    while (!q.empty())
    {
        // std::cout << "BFS" << std::endl;
        int current = q.front();
        q.pop();
        largestComponent.insert(current);
        for (int neighbor : graphForCore[current])
        {
            if (visited.find(neighbor) == visited.end() && core[neighbor] >= k)
            {
                q.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }
    return largestComponent;
}

int main(int argc, char *argv[])
{
    // 检查是否有输入文件路径
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <path_to_dataset>" << argv[1] << " if use opt_index cons?" << argv[2] << " if use serialize_index?" << std::endl;
        return 1;
    }
    // 获取命令行输入的文件路径
    std::string dataset_path = argv[1];
    std::string indexFile = dataset_path + ".index"; // 索引文件路径
    std::string opt_flag = argv[2];
    std::string index_flag = argv[3];

    Graph graph1(dataset_path);

    std::cout << "图加载完成，开始构造标签集合以及查询标签" << std::endl;

    std::unordered_set<std::string> lables = graph1.checkLabelForm();
    std::unordered_set<std::string> queryLabels;

    // LCS optLcs1 = LCS(lables, graph1, 1);
    // optLcs1.printLcsIndex();
    // return 0;

    if (lables.size() == 8)
    {
        if (lables.find("0") != lables.end())
        {
            queryLabels = {"0", "2", "4", "6"};
        }
        else
        {
            queryLabels = {"1", "3", "5", "7"};
        }
    }
    else if (lables.size() == 5)
    {
        queryLabels = {"0", "2"};
    }
    else if (lables.size() == 20)
    {
        queryLabels = {"0", "2", "4", "6", "8", "10", "12", "14", "16", "18"};
    }

    std::cout << "选定标签集合和查询标签集合" << std::endl;

    int count = 3;

    // 随机数引擎
    std::random_device rd;
    std::mt19937 gen(rd());

    // 开始得到连通图
    std::cout << "查看当前数据集是否为连通图" << std::endl;

    // BFS
    std::unordered_set<int> nodes = graph1.findLargestConnectedComponent();
    std::unordered_set<int> vertexs = graph1.getNodes();
    std::vector<int> nodesToRemove;
    int sum = nodes.size();
    if (sum != vertexs.size())
    {
        std::cout << "Graph is not connected." << std::endl;
        // 找出所有不在最大连通分量中的节点
        for (const auto &node : vertexs)
        {
            if (nodes.find(node) == nodes.end())
            {
                nodesToRemove.push_back(node);
            }
        }

        // 删除这些节点
        for (int node : nodesToRemove)
        {
            graph1.removeNode(node);
        }
    }
    else
    {
        std::cout << "Graph is connected." << std::endl;
    }

    std::cout << "make the graph connected." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    Hierarchy hierarchy(graph1);
    int Kmax = hierarchy.printVerticesByCoreness();

    std::cout << "successedfully build the k-core hierarchy." << std::endl;

    auto lcsStart = std::chrono::high_resolution_clock::now();

    std::cout << "the k-core hierarchy built in " << std::chrono::duration_cast<std::chrono::milliseconds>(lcsStart - start).count() << " ms" << std::endl;

    LCS optLcs; // 创建LCS对象

    // 检查索引文件是否已经存在
    if (index_flag == "true")
    {
        if (std::filesystem::exists(indexFile))
        {
            std::cout << "Loading index from file..." << std::endl;
            if (!deserializeIndex(optLcs, indexFile))
            {
                std::cerr << "Failed to load index." << std::endl;
                return 1;
            }
        }
        else
        {
            std::cout << "Index file not found. Building index..." << std::endl;
            // 如果索引文件不存在，构建索引
            if (opt_flag == "true")
            {
                optLcs = LCS(lables, graph1, 1);
            }
            else if (opt_flag == "false")
            {
                optLcs = LCS(lables, graph1, 0);
            }
            std::cout << "Saving index to file..." << std::endl;
            serializeIndex(optLcs, indexFile); // 保存索引到文件
        }
    }
    else
    {
        if (opt_flag == "true")
        {
            optLcs = LCS(lables, graph1, 1);
        }
        else if (opt_flag == "false")
        {
            optLcs = LCS(lables, graph1, 0);
        }
    }

    std::cout << "the LCS index has been built." << std::endl;

    auto mid = std::chrono::high_resolution_clock::now();

    std::cout << "the LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;

    // TreeIndex index = TreeIndex(graph1, lables);

    // std::cout << "the K-Shell index has been built." << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    double indexCost[5] = {0, 0, 0, 0, 0};
    double onlineCost[5] = {0, 0, 0, 0, 0};
    double vaildQuery[5] = {0, 0, 0, 0, 0};
    double baseline[5] = {0, 0, 0, 0, 0};

    double anti_indexCost[5] = {0, 0, 0, 0, 0};
    double anti_onlineCost[5] = {0, 0, 0, 0, 0};
    double anti_vaildQuery[5] = {0, 0, 0, 0, 0};

    // case study
    // std::unordered_set<int> solution1 = hierarchy.query(query, k);

    for (int i = 0; i < 5; i++)
    {
        std::vector<int> validNumbers = hierarchy.vaildCandicate(i, Kmax);

        static std::mt19937 gen(std::time(0));

        int k = 0.1 * (i * 2 + 1) * Kmax;
        if (k == 1)
            k++;
        // note:此处是从层次结构的k挑选顶点 因此在边标签约束下可能已经达不到此k 以此模拟输入不符合查询的情况
        std::uniform_int_distribution<> dis(0, validNumbers.size() - 1);

        int vaild = 0;

        int invaild = 0;

        int anti_vaild = 0;

        for (int j = 0; j < 10; j++)
        {

            int query = validNumbers[dis(gen)];

            std::vector<int> queryNodes;

            queryNodes.push_back(query);

            auto begin = std::chrono::high_resolution_clock::now();

            Graph temp1(graph1);

            std::unordered_set<int> solution1 = hierarchy.bl(query, k, temp1, queryLabels);

            auto query1 = std::chrono::high_resolution_clock::now();

            // std::unordered_set<int> solution2 = CommunitySearch::online(graph1, queryNodes, k);

            auto query2 = std::chrono::high_resolution_clock::now();

            // std::unordered_set<int> solution3 = index.findKCoreSubgraph(queryNodes, k);

            auto query3 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution4 = hierarchy.queryWithLcs(graph1, query, k, optLcs, queryLabels);

            auto query4 = std::chrono::high_resolution_clock::now();

            Graph temp(graph1);

            auto queryt = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution5 = CommunitySearch::onlineWithEdge(temp, queryNodes, k, queryLabels);

            auto query5 = std::chrono::high_resolution_clock::now();

            // std::cout<<"arrive";

            // Graph temp1(graph1);

            // std::cout<<"arrive";

            // temp1.filterEdgeLableAndRecover(queryLabels);

            // std::unordered_set<int> solution6 = computeCore(temp1,*queryNodes.begin(),k);

            // std::cout << "the size of solution1 is " << solution1.size() << std::endl;
            // std::cout << "the size of solution2 is " << solution2.size() << std::endl;
            // std::cout << "the size of solution3 is " << solution3.size() << std::endl;
            // std::cout << "the size of solution4 is " << solution4.size() << std::endl;
            // std::cout << "the size of solution5 is " << solution5.size() << std::endl;
            if (solution4.size() == solution5.size() && solution1.size() == solution4.size() && solution4.size() >= k)
            {
                indexCost[i] += std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count();
                // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count() << std::endl;
                onlineCost[i] += std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count();
                // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count() << std::endl;
                vaildQuery[i] += solution4.size();
                vaild++;
                baseline[i]+=std::chrono::duration_cast<std::chrono::milliseconds>(query1 - begin).count();
            }
            else
            {
                anti_indexCost[i] += std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count();
                // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count() << std::endl;
                anti_onlineCost[i] += std::chrono::duration_cast<std::chrono::milliseconds>(query5 - queryt).count();
                // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count() << std::endl;
                anti_vaildQuery[i] += solution4.size();
                invaild++;
                if (!optLcs.isExist(query, queryLabels, k) && solution5.size() == 0)
                {
                    anti_vaild++;
                }
                else
                {
                    // std::cout<<"the index query solution's size is "<<solution4.size()<<std::endl;
                    // std::cout<<"the online query solution's size is "<<solution5.size()<<std::endl;
                }
                // continue;
            }
        }
        if (vaild != 0)
        {
            indexCost[i] /= vaild;
            onlineCost[i] /= vaild;
            vaildQuery[i] /= vaild;
            baseline[i] /= vaild;

            std::cout << "k为 " << k << " 时,基于索引的平均查询时间为" << indexCost[i] << std::endl;
            std::cout << "k为 " << k << " 时,在线搜索的平均查询时间为" << onlineCost[i] << std::endl;
            std::cout << "k为 " << k << " 时,层次结构的平均查询时间为" << baseline[i] << std::endl;
            std::cout << "k为 " << k << " 时,平均查询结果数量为" << vaildQuery[i] << std::endl;
            std::cout << "k为 " << k << " 时,查询成功的次数为" << vaild << std::endl;
        }
        if (invaild != 0)
        {
            anti_indexCost[i] /= invaild;
            anti_onlineCost[i] /= invaild;
            anti_vaildQuery[i] /= invaild;

            std::cout << "k为 " << k << " 时,基于索引的平均查询时间为" << anti_indexCost[i] << std::endl;
            std::cout << "k为 " << k << " 时,在线搜索的平均查询时间为" << anti_onlineCost[i] << std::endl;
            std::cout << "k为 " << k << " 时,平均查询结果数量为" << anti_vaildQuery[i] << std::endl;
            std::cout << "k为 " << k << " 时,查询失败的次数为" << invaild << std::endl;
            std::cout << "k为 " << k << " 时,query不规范次数为" << anti_vaild << std::endl;
        }
        std::cout << std::endl;
    }

    // // 判断是否相同
    // int flag = 1;

    // if (flag)
    // {
    //     std::cout << "PASS!!!!!!!!!!!!!!!!!" << std::endl;
    //     // std::cout << "the normal LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(optLcsStart - lcsStart).count() << " ms" << std::endl;
    //     std::cout << "the optLCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;
    //     // std::cout << "the k-shell struct built in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count() << " ms" << std::endl<< std::endl;
    //     std::cout << "the k-core hierarchy query without label in " << std::chrono::duration_cast<std::chrono::milliseconds>(query1 - begin).count() << " ms" << std::endl;
    //     std::cout << "the online search query without label in " << std::chrono::duration_cast<std::chrono::milliseconds>(query2 - query1).count() << " ms" << std::endl;
    //     // std::cout << "the k-shell struct query in " << std::chrono::duration_cast<std::chrono::milliseconds>(query3 - query2).count() << " ms" << std::endl<< std::endl;;
    //     std::cout << "the k-core hierarchy egde query in " << std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count() << " ms" << std::endl;
    //     std::cout << "the online search edge query in " << std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count() << " ms" << std::endl;
    // }

    /************************************************************/

    // TreeIndex index = TreeIndex(graph1, lables);

    // auto mid = std::chrono::high_resolution_clock::now();
    //  计算时间差
    // auto midduration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start);
    //  输出索引建立运行时间
    // std::cout << "索引建立算法运行时间: " << midduration.count() << " 毫秒" << std::endl;

    // index.printQueryNodesCoreIndex(queryNodes);
    // std::unordered_set<int> queryNodesSet(queryNodes.begin(), queryNodes.end());

    // if (GraphUtility::isConnected(queryNodesSet, graph))
    // {
    //     std::cout << "Query nodes are connected." << std::endl;
    // }
    // else
    // {
    //     std::cout << "Query nodes are not connected." << std::endl;
    // }

    // auto newStart = std::chrono::high_resolution_clock::now();

    // solution1 = index.findKCoreSubgraph(queryNodes, k);

    // auto end = std::chrono::high_resolution_clock::now();

    // // 计算时间差
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - newStart);

    // // 输出程序运行时间
    // std::cout << "基于索引重建连通分量算法运行时间: " << duration.count() << " 毫秒" << std::endl;

    // std::cout << solution1.size() << std::endl;
    // if (GraphUtility::isConnectedWithLabels(solution1, graph1, lables))
    // {
    //     std::cout << "Solution 1 is valid." << std::endl;
    // }
    // else
    // {
    //     std::cout << "Solution 1 is invalid." << std::endl;
    // }
    // std::unordered_map<int, int> core = CoreGroup::coreGroupsAlgorithmWithEdge(graph1, lables, true);
    // for (int i : solution1)
    // {
    //     if (core[i] < 4)
    //     {
    //         std::cout << "error" << std::endl;
    //         break;
    //     }
    // }
    // if (GraphUtility::testSolution(solution1, queryNodes, index, graph, k))
    // {
    //     std::cout << "Solution 1 is valid." << std::endl;
    // }
    // else
    // {
    //     std::cout << "Solution 1 is invalid." << std::endl;
    // }

    // // 计算在线搜索算法
    // start = std::chrono::high_resolution_clock::now();
    // Graph graph2("./Dataset/amazon-meta-graph-label.txt");
    // solution2 = CommunitySearch::onlineWithEdge(graph2, queryNodes, k, lables);

    // end = std::chrono::high_resolution_clock::now();

    // 计算时间差
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // 输出程序运行时间

    // std::cout << "在线搜索算法运行时间: " << duration.count() << " 毫秒" << std::endl;

    // std::cout << solution2.size() << std::endl;

    // if (GraphUtility::isConnectedWithLabels(solution2, graph2, lables))
    // {
    //     std::cout << "Solution 2 is valid." << std::endl;
    // }
    // else
    // {
    //     std::cout << "Solution 2 is invalid." << std::endl;
    // }

    // core = CoreGroup::coreGroupsAlgorithmWithEdge(graph2, lables, true);
    // for (int i : solution2)
    // {
    //     if (core[i] < k)
    //     {
    //         std::cout << "error" << std::endl;
    //         break;
    //     }
    // }

    // for(auto node:solution2){
    //     std::cout<<node<<std::endl;
    // }
    // if (GraphUtility::testSolution(solution2, queryNodes, index, graph, k))
    // {
    //     std::cout << "Solution 2 is valid." << std::endl;
    // }
    // else
    // {
    //     std::cout << "Solution 2 is invalid." << std::endl;
    // }

    // std::unordered_set<int> difference;
    // for (int node : solution2)
    // {
    //     if (solution1.find(node) == solution1.end())
    //     {
    //         difference.insert(node);
    //     }
    // }
    // std::unordered_map<int, int> coreIndexMap;
    // for (int node : difference)
    // {
    //     coreIndexMap[node] = index.getCoreIndex(node);
    // }
    // std::cout << "Core index for elements in the difference set:" << std::endl;
    // for (const auto &pair : coreIndexMap)
    // {
    //     std::cout << "Node: " << pair.first << ", Core index: " << pair.second << std::endl;
    // }

    // // 1. 找到这几个节点所属的连通分量的 ID
    // std::unordered_map<int, int> componentIds;
    // for (int node : difference)
    // {
    //     int componentId = index.nodeToComponentId[node];
    //     componentIds[node] = componentId;
    // }

    // // 2. 打印这几个节点所属的连通分量的父子索引
    // for (const auto &pair : componentIds)
    // {
    //     int node = pair.first;
    //     int componentId = pair.second;

    //     std::cout << "Node " << node << " belongs to component " << componentId << std::endl;

    //     int parentComponentId = index.getParentComponent(componentId);
    //     std::cout << "Parent component of component " << componentId << " is " << parentComponentId << std::endl;

    //     std::unordered_set<int> childrenComponents = index.getConnectedComponentChildren(componentId);
    //     std::cout << "Children components of component " << componentId << " are: ";
    //     for (int child : childrenComponents)
    //     {
    //         std::cout << child << " ";
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}
class GraphUtility
{
public:
    // 测试solution中的节点的最小k是否大于等于查询节点组的k，并检查连通性
    static bool testSolution(const std::unordered_set<int> &solution, const std::vector<int> &queryNodes, TreeIndex &index, Graph &graph, int k)
    {
        if (solution.empty())
            return false;

        // 获取solution中的最小核心级别
        int minCoreIndexSolution = INT_MAX;
        for (int node : solution)
        {
            int coreIndex = index.getCoreIndex(node);
            minCoreIndexSolution = std::min(minCoreIndexSolution, coreIndex);
        }
        std::cout << "Solution's min k is :" << minCoreIndexSolution << std::endl;

        // 获取查询节点的最小核心级别
        int minCoreIndexQuery = INT_MAX;
        for (int node : queryNodes)
        {
            int coreIndex = index.getCoreIndex(node);
            minCoreIndexQuery = std::min(minCoreIndexQuery, coreIndex);
        }

        // 判断最小核心级别
        if (k > minCoreIndexSolution)
        {
            std::cout << "Solution's min k is less than the given k." << std::endl;
            return false;
        }

        // 检查连通性
        std::unordered_set<int> solutionSet(solution.begin(), solution.end());
        return isConnected(solutionSet, graph);
    }

    // 使用BFS检查一组节点是否连通
    static bool isConnected(const std::unordered_set<int> &nodes, Graph &graph)
    {
        if (nodes.empty())
            return false;
        std::unordered_set<int> visited;
        std::unordered_set<int> temp = nodes;
        std::queue<int> queue;

        // 从任意一个节点开始
        int startNode = *nodes.begin();
        queue.push(startNode);
        visited.insert(startNode);

        // 进行BFS
        while (!queue.empty())
        {
            int current = queue.front();
            queue.pop();
            // 从temp中删除已经访问过的节点
            temp.erase(current);
            for (int neighbor : graph.getNeighbors(current))
            {
                if (nodes.find(neighbor) != nodes.end() && visited.find(neighbor) == visited.end())
                {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        // 检查是否所有节点都被访问过
        return temp.size() == 0;
    }

    static bool isConnectedWithLabels(const std::unordered_set<int> &nodes, Graph &graph, std::unordered_set<std::string> lables)
    {
        if (nodes.empty())
            return false;
        std::unordered_set<int> visited;
        std::unordered_set<int> temp = nodes;
        std::queue<int> queue;

        // 从任意一个节点开始
        int startNode = *nodes.begin();
        queue.push(startNode);
        visited.insert(startNode);

        // 进行BFS
        while (!queue.empty())
        {
            int current = queue.front();
            queue.pop();
            // 从temp中删除已经访问过的节点
            temp.erase(current);
            for (int neighbor : graph.getNeighbors(current))
            {
                for (std::string lable : lables)
                {
                    if (graph.getEdgeLable(current, neighbor).find(lable) != graph.getEdgeLable(current, neighbor).end())
                    {
                        if (nodes.find(neighbor) != nodes.end() && visited.find(neighbor) == visited.end())
                        {
                            visited.insert(neighbor);
                            queue.push(neighbor);
                        }
                        break;
                    }
                }
            }
        }

        // 检查是否所有节点都被访问过
        return temp.size() == 0;
    }
};