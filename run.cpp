#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <string.h>
#include <unordered_set>
#include <random>
#include <ctime>
#include "Graph.h"
#include "./Algorithm/CoreGroup.h"
#include "./BuildIndex/TreeIndex.h"
#include "./BuildIndex/hierarchy.h"
#include "./BuildIndex/LCS.h"
#include "./Algorithm/Online.h"
#include "./utils/Mapping.h"

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

int main(int argc, char *argv[])
{
    // 创建图,并且初始化数据,输入图必须连通,且无向

    Graph graph1("/home/asc23/lcy/graph/data/semanticscholar/case_study/sem_network8.txt");

    std::unordered_set<std::string> lables = {"0", "1", "2", "3", "4", "5", "6", "7"};
    std::unordered_set<std::string> queryLabels = {"0", "2","5"};

    // BFS
    std::vector<int> queryNodes;

    queryNodes.push_back(249);
    for (int i = 5; i <= 10; i++)
    {
        Graph tempgraph(graph1);
        std::unordered_set<int> solution4 = CommunitySearch::onlineWithEdge(tempgraph, queryNodes, i, queryLabels);
        std::cout<<"the solution is "<<solution4.size()<<std::endl;
        // 寻找
        std::unordered_set<int> solution5;
        std::queue<int> q1;
        std::unordered_set<int> visited1;
        q1.push(249);
        visited1.insert(249);
        for (int j = 0; j < 3; j++)
        {
            std::queue<int> temp = q1;
            while (!temp.empty())
            {
                int w = temp.front();
                solution5.insert(w);
                temp.pop();
                q1.pop();
                for (int v : graph1.getNeighbors(w))
                {
                    if (visited1.find(v) == visited1.end() && solution4.find(v) != solution4.end())
                    {
                        visited1.insert(v);
                        q1.push(v);
                    }
                }
            }
        }
        for (int j : solution5)
        {
            std::cout << j << " ";
        }

        std::cout << std::endl;
        std::cout << "query with label constraints' answer as above ; and the sum is " << solution5.size() << std::endl;
    }

    // auto start = std::chrono::high_resolution_clock::now();

    // Hierarchy hierarchy(graph1);
    // int Kmax = hierarchy.printVerticesByCoreness();

    // std::cout << "successedfully build the k-core hierarchy." << std::endl;

    // auto lcsStart = std::chrono::high_resolution_clock::now();

    // LCS optLcs(lables, graph1);

    // std::cout << "the LCS index has been built." << std::endl;

    // auto mid = std::chrono::high_resolution_clock::now();

    // // TreeIndex index = TreeIndex(graph1, lables);

    // std::cout << "the K-Shell index has been built." << std::endl;

    // auto end = std::chrono::high_resolution_clock::now();

    // double indexCost[5];
    // double onlineCost[5];

    // // case study
    // // std::unordered_set<int> solution1 = hierarchy.query(query, k);

    // for (int i = 0; i < 5; i++)
    // {
    //     std::vector<int> validNumbers = hierarchy.vaildCandicate(i, Kmax);

    //     static std::mt19937 gen(std::time(0));

    //     int k = 0.2 * (i + 1) * Kmax;

    //     std::uniform_int_distribution<> dis(0, validNumbers.size() - 1);

    //     for (int j = 0; j < 1000; j++)
    //     {
    //         int query = validNumbers[dis(gen)];

    //         std::vector<int> queryNodes;

    //         queryNodes.push_back(query);

    //         auto begin = std::chrono::high_resolution_clock::now();

    //         // std::unordered_set<int> solution1 = hierarchy.query(query, k);

    //         auto query1 = std::chrono::high_resolution_clock::now();

    //         // std::unordered_set<int> solution2 = CommunitySearch::online(graph1, queryNodes, k);

    //         auto query2 = std::chrono::high_resolution_clock::now();

    //         // std::unordered_set<int> solution3 = index.findKCoreSubgraph(queryNodes, k);

    //         auto query3 = std::chrono::high_resolution_clock::now();

    //         std::unordered_set<int> solution4 = hierarchy.queryWithLcs(graph1, query, k, optLcs, queryLabels);

    //         auto query4 = std::chrono::high_resolution_clock::now();

    //         std::unordered_set<int> solution5 = CommunitySearch::onlineWithEdge(graph1, queryNodes, k, queryLabels);

    //         auto query5 = std::chrono::high_resolution_clock::now();

    //         // std::cout << "the size of solution1 is " << solution1.size() << std::endl;
    //         // std::cout << "the size of solution2 is " << solution2.size() << std::endl;
    //         // std::cout << "the size of solution3 is " << solution3.size() << std::endl;
    //         // std::cout << "the size of solution4 is " << solution4.size() << std::endl;
    //         // std::cout << "the size of solution5 is " << solution5.size() << std::endl;
    //         if (solution4.size() == solution5.size())
    //         {
    //             indexCost[i - 1] += std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count();
    //             // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count() << std::endl;
    //             onlineCost[i - 1] += std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count();
    //             // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(query5 - query4).count() << std::endl;
    //         }
    //         else
    //         {
    //             continue;
    //         }
    //     }
    //     indexCost[i - 1] /= 1000;
    //     onlineCost[i - 1] /= 1000;
    //     std::cout << "k为 " << k << " 时,基于索引的平均查询时间为" << indexCost[i - 1] << std::endl;
    //     std::cout << "k为 " << k << " 时,在线搜索的平均查询时间为" << onlineCost[i - 1] << std::endl;
    // }

    // // 判断是否相同
    // int flag = 1;

    // if (flag)
    // {
    //     std::cout << "PASS!!!!!!!!!!!!!!!!!" << std::endl;
    // std::cout << "the k-core hierarchy built in " << std::chrono::duration_cast<std::chrono::milliseconds>(lcsStart - start).count() << " ms" << std::endl;
    //     // std::cout << "the normal LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(optLcsStart - lcsStart).count() << " ms" << std::endl;
    // std::cout << "the optLCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;
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
