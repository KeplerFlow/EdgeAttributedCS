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

// Helper function to generate combinations
void generateCombinations(const std::vector<std::string> &labels, int combinationSize, int startIndex, std::vector<std::string> &currentCombination, std::vector<std::unordered_set<std::string>> &allCombinations)
{
    if (currentCombination.size() == combinationSize)
    {
        allCombinations.push_back(std::unordered_set<std::string>(currentCombination.begin(), currentCombination.end()));
        return;
    }

    for (int i = startIndex; i < labels.size(); ++i)
    {
        currentCombination.push_back(labels[i]);
        generateCombinations(labels, combinationSize, i + 1, currentCombination, allCombinations);
        currentCombination.pop_back();
    }
}

// Function to get all combinations of labels from size 1 to maxCombinationSize
std::vector<std::unordered_set<std::string>> getAllLabelCombinations(const std::unordered_set<std::string> &labelSet, int maxCombinationSize)
{
    std::vector<std::string> labels(labelSet.begin(), labelSet.end());
    std::vector<std::unordered_set<std::string>> allCombinations;

    for (int i = 1; i <= maxCombinationSize; ++i)
    {
        std::vector<std::string> currentCombination;
        generateCombinations(labels, i, 0, currentCombination, allCombinations);
    }

    return allCombinations;
}

std::unordered_set<int> queryWithLcs(Graph &graph, int k, LCS &lcs, std::unordered_set<std::string> &labels, std::vector<std::unordered_set<int>> &coreness, int kmax)
{
    std::unordered_set<int> solution;
    for (int i = k; i <= kmax; i++)
    {
        if (!coreness[i].empty())
        {
            for (int node : coreness[i])
            {
                if (lcs.isExist(node, labels, k))
                {
                    solution.insert(node);
                }
            }
        }
    }
    return solution;
}

std::unordered_set<int> queryWithMuptIndex(Graph &graph, int k, LCS &lcs, std::unordered_set<std::string> &labels)
{
    std::unordered_set<int> solution;
    std::unordered_map<int,std::unordered_set<int>> core = lcs.mupt_index_search(labels);
    for(auto &pair:core){
        if(pair.first>=k){
            solution.insert(pair.second.begin(),pair.second.end());
        }
    }
    return solution;
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

    Graph graph1(dataset_path,2336696);

    std::cout << "图加载完成，开始构造标签集合以及查询标签" << std::endl;

    std::unordered_set<std::string> lables = graph1.checkLabelForm();
    std::vector<std::unordered_set<std::string>> allCombinations;

    int maxCombinationSize = lables.size();
    allCombinations = getAllLabelCombinations(lables, maxCombinationSize);

    std::cout << "选定标签集合和查询标签集合" << std::endl;

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

    std::unordered_map<int, int> coreIndex = CoreGroup::coreGroupsAlgorithm(graph1);
    // std::unordered_map<int, std::unordered_set<int>> coreness;
    std::vector<std::unordered_set<int>> core;
    int kmax = 0;

    for (auto &pair : coreIndex)
    {
        // std::cout << "the vertex is " << pair.first << " the coreness is " << pair.second << std::endl;
        // coreness[pair.second].insert(pair.first);
        if (pair.second > kmax)
        {
            kmax = pair.second;
            core.resize(kmax + 1);
        }
        core[pair.second].insert(pair.first);
    }

    //auto lcsStart = std::chrono::high_resolution_clock::now();

    // std::cout << "the k-core hierarchy built in " << std::chrono::duration_cast<std::chrono::milliseconds>(lcsStart - start).count() << " ms" << std::endl;

    LCS optLcs(lables, graph1); // 创建LCS对象

    //auto mid = std::chrono::high_resolution_clock::now();

    //std::cout << "the basic LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;
  
    LCS lcs(lables, graph1,1);
    std::unordered_set<std::string> queryLabels = {"0", "2", "4", "6"};
    long long indexCost=0,onlineCost=0;
    int vaild=0;
    {
        for(int i=0;i<10;i++)
        {
            int k = 0.6 * kmax;

            if (k == 1)
                k++;
            else if (k == kmax)
                k--;

            auto query1 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution4 = queryWithMuptIndex(graph1, k, optLcs, queryLabels);

            auto query2 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution5 = queryWithLcs(graph1, k, lcs, queryLabels, core, kmax);

            auto query3 = std::chrono::high_resolution_clock::now();

            if(solution4.size()==solution5.size()&&solution4.size()>1)
            {
                indexCost += std::chrono::duration_cast<std::chrono::nanoseconds>(query2 - query1).count();
                onlineCost += std::chrono::duration_cast<std::chrono::nanoseconds>(query3 - query2).count();
                vaild++;
            }
        }
        std::cout <<  "基础索引的查询时间为" << indexCost/ vaild<<" ns" <<std::endl;
        std::cout << "优化索引的查询时间为" << onlineCost / vaild<< " ns"<<std::endl;
        std::cout << std::endl;

    }

}