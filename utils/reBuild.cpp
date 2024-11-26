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
    Graph graph2(dataset_path,100,120,"7");

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

    nodes = graph2.findLargestConnectedComponent();
    vertexs = graph2.getNodes();
    nodesToRemove;
    sum = nodes.size();
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
            graph2.removeNode(node);
        }
    }
    else
    {
        std::cout << "Graph is connected." << std::endl;
    }

    std::cout << "make the graph connected." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

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
    //对比
    LCS comLcs = LCS(lables, graph2, 1);

    //需要增量维护
    optLcs.incIndex(graph1,100,120,"7");
    //comLcs.incIndex(graph1,100,120,"7");

    std::unordered_map<int, std::unordered_map<int, std::vector<std::unordered_set<std::string>>>> lhs = optLcs.getLcsIndex();
    std::unordered_map<int, std::unordered_map<int, std::vector<std::unordered_set<std::string>>>> rhs = comLcs.getLcsIndex();

    std::cout<<"begin to compare the index value"<<std::endl;

    for (const auto& outer_pair : lhs) {
        const int outer_key = outer_pair.first;
        const auto& lhs_inner_map = outer_pair.second;

        // Check if the key exists in the second map
        auto it = rhs.find(outer_key);
        if (it == rhs.end()){
            std::cerr<<"最外层顶点不同";
            return 0;
        }

        const auto& rhs_inner_map = it->second;

        // Compare inner maps
        if (lhs_inner_map.size() != rhs_inner_map.size()) 
        {
            std::cerr<<"次外层k数目不同";
            return 0;
        }

        for (const auto& inner_pair : lhs_inner_map) {
            const int inner_key = inner_pair.first;
            const auto& lhs_vector = inner_pair.second;

            // Check if the key exists in the second inner map
            auto inner_it = rhs_inner_map.find(inner_key);
            if (inner_it == rhs_inner_map.end())
            {
                std::cerr<<"次外层k值不同";
                return 0;
            }

            const auto& rhs_vector = inner_it->second;

            // Compare vectors
            if (lhs_vector.size() != rhs_vector.size()){
                std::cerr<<"内层标签集合数目不同";
                return 0;
            }

            for (size_t i = 0; i < lhs_vector.size(); ++i) {
                // Compare unordered sets
                if (lhs_vector[i] != rhs_vector[i]) {
                    std::cerr<<"内层标签不同";
                    return 0;
                }
            }
        }
    }
    std::cerr<<"PASS!!"<<std::endl;

}

