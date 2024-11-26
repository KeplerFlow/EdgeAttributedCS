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
    /**** */
    // std::unordered_map<int,std::unordered_set<int>> cores = lcs.mupt_index_search(labels);
    // for(auto &pair:cores){
    //     if(pair.first>=k){
    //         solution.insert(pair.second.begin(),pair.second.end());
    //     }
    // }
    /**** */
    return lcs.mupt_index_search_map(labels,k);
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
    std::string muptindexFile = dataset_path + "_index.bin"; // 动态生成文件名
    std::string opt_flag = argv[2];
    std::string index_flag = argv[3];

    Graph graph1(dataset_path);

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

    auto start = std::chrono::high_resolution_clock::now();

    // Hierarchy hierarchy(graph1);
    // int Kmax = hierarchy.printVerticesByCoreness();

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

    std::cout << "the K max is " << kmax << std::endl;

    std::cout << "successedfully build the k-core ." << std::endl;

    auto lcsStart = std::chrono::high_resolution_clock::now();

    // std::cout << "the k-core hierarchy built in " << std::chrono::duration_cast<std::chrono::milliseconds>(lcsStart - start).count() << " ms" << std::endl;

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

    auto mid = std::chrono::high_resolution_clock::now();

    std::cout << "the LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;
    size_t memoryUsage = optLcs.calculateMemoryUsage();
    std::cout << "Memory usage of LcsIndex: " << memoryUsage / (1024.0 * 1024.0) << " MB" << std::endl;

    auto basic = std::chrono::high_resolution_clock::now();

    LCS basicLcs; // 创建basic LCS对象
    if (index_flag == "true") {
        if (std::filesystem::exists(muptindexFile)) {
            std::cout << "Loading basic index from file: " << muptindexFile << std::endl;
            std::ifstream ifs(muptindexFile, std::ios::binary);
            if (!ifs)
            {
                std::cerr << "Error opening file for reading: " << muptindexFile << std::endl;
                return false;
            }
            basicLcs.deserialize_mupt_index(ifs);
            basicLcs.enumLabels(lables);
            ifs.close();
        } else {
            std::cout << "Index file not found. Building index for dataset: " << dataset_path << std::endl;

            // 如果索引文件不存在，构建索引
            basicLcs = LCS(lables, graph1);

            std::cout << "Saving index to file: " << muptindexFile << std::endl;
            std::ofstream ofs(muptindexFile, std::ios::binary);
            if (!ofs)
            {
                std::cerr << "Error opening file for writing: " << muptindexFile << std::endl;
                return 0;
            }
            // 这里假设LCS类可以被序列化（可以自己实现LCS的序列化函数）
            basicLcs.serialize_mupt_index(ofs);
            ofs.close();
        }
    }else{
        basicLcs = LCS(lables, graph1);
    }

    // std::cout << "the basic index has been built." << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    //std::cout << "the basic LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - basic).count() << " ms" << std::endl;

    bool query_label = false;
    basicLcs.transIndex();
    if(false)
    {
        std::cout<<"进行默认值实验"<<std::endl;
        std::unordered_set<std::string> queryLabels;
        if (lables.size() == 8)
        {
            if (lables.find("0") != lables.end())
            {
                queryLabels = {"0", "2", "4", "6"};
            }
            else if(lables.find("17") != lables.end())
            {
                queryLabels = {"5", "7", "14"};
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
        else if(lables.size() ==10)
        {
            if(lables.find("15") != lables.end())
            {
                queryLabels = {"0", "4", "10", "14", "5"}; 
            }
            else
            {
                queryLabels = {"0", "2", "4", "6", "8"};
            }
        }else if(lables.size()==12){
            queryLabels = {"0", "2", "4", "6", "8", "10"};
        }else if(lables.size()==11){
            queryLabels = {"2", "3","6","7","9"};
        }
        int k = 0.6*kmax;
        double indexCost = 0, basicCost = 0;
        for (int j = 0; j < 100; j++) {
            auto query1 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution4 = queryWithLcs(graph1, k, optLcs, queryLabels, core, kmax);

            auto query2 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution5 = queryWithMuptIndex(graph1, k, basicLcs, queryLabels);

            auto query3 = std::chrono::high_resolution_clock::now();

            if (solution4.size() == solution5.size()) {
                indexCost += std::chrono::duration_cast<std::chrono::nanoseconds>(query2 - query1).count();
                basicCost += std::chrono::duration_cast<std::chrono::nanoseconds>(query3 - query2).count();
            } else {
                std::cout << "The size of LCS index is " << solution4.size() 
                        << ", the size of basic index is " << solution5.size() << std::endl;
            }
        }

        std::cout << "Average query time for index: " << indexCost / 100 << " ns" << std::endl;
        std::cout << "Average query time for basic index: " << basicCost / 100 << " ns" << std::endl;
        return 0;
    }

    // case study
    // std::unordered_set<int> solution1 = hierarchy.query(query, k);
    if (query_label == true)
    {
        std::cout << "进行标签调整实验" << std::endl;
        int exp_1 = 0.2 * lables.size();
        int exp_2 = 0.4 * lables.size();
        int exp_3 = 0.6 * lables.size();
        int exp_4 = 0.8 * lables.size();
        double indexCost[exp_4 + 1] = {0};
        double onlineCost[exp_4 + 1] = {0};
        int vaild[exp_4 + 1] = {0};

        for (auto &queryLabels : allCombinations)
        {
            if (queryLabels.size() != exp_1 && queryLabels.size() != exp_2 && queryLabels.size() != exp_3 && queryLabels.size() != exp_4)
            {
                continue;
            }
            int k = 0.6 * kmax;

            if (k == 1)
                k++;
            else if (k == kmax)
                k--;

            auto query1 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution4 = queryWithLcs(graph1, k, optLcs, queryLabels, core, kmax);

            auto query2 = std::chrono::high_resolution_clock::now();

            //Graph temp(graph1);
            //std::unordered_set<int> solution5 = CommunitySearch::onlineWithEdge(temp, k, queryLabels);
            
            std::unordered_set<int> solution5 = queryWithMuptIndex(graph1, k, basicLcs, queryLabels);

            auto query3 = std::chrono::high_resolution_clock::now();

            if(solution4.size()==solution5.size())
            {
                indexCost[queryLabels.size()] += std::chrono::duration_cast<std::chrono::milliseconds>(query2 - query1).count();
                onlineCost[queryLabels.size()] += std::chrono::duration_cast<std::chrono::milliseconds>(query3 - query2).count();
                vaild[queryLabels.size()]++;
            }
        }
        std::cout << "标签大小为 " << exp_1 << " 时,基于索引的查询时间为" << indexCost[exp_1] / vaild[exp_1] << std::endl;
        std::cout << "标签大小为 " << exp_1 << " 时,对照算法的查询时间为" << onlineCost[exp_1] / vaild[exp_1] << std::endl;
        std::cout << std::endl;
        std::cout << "标签大小为 " << exp_2 << " 时,基于索引的查询时间为" << indexCost[exp_2] / vaild[exp_2] << std::endl;
        std::cout << "标签大小为 " << exp_2 << " 时,对照算法的查询时间为" << onlineCost[exp_2] / vaild[exp_2] << std::endl;
        std::cout << std::endl;
        std::cout << "标签大小为 " << exp_3 << " 时,基于索引的查询时间为" << indexCost[exp_3] / vaild[exp_3] << std::endl;
        std::cout << "标签大小为 " << exp_3 << " 时,对照算法的查询时间为" << onlineCost[exp_3] / vaild[exp_3] << std::endl;
        std::cout << std::endl;
        std::cout << "标签大小为 " << exp_4 << " 时,基于索引的查询时间为" << indexCost[exp_4] / vaild[exp_4] << std::endl;
        std::cout << "标签大小为 " << exp_4 << " 时,对照算法的查询时间为" << onlineCost[exp_4] / vaild[exp_4] << std::endl;
        std::cout << std::endl;

    }
    else
    {
        std::cout << "进行k值调整实验" << std::endl;
        std::unordered_set<std::string> queryLabels;
        double indexCost = 0;
        double onlineCost = 0;
        if (lables.size() == 8)
        {
            if (lables.find("0") != lables.end())
            {
                queryLabels = {"0", "2", "4", "6"};
            }
            else if(lables.find("17") != lables.end())
            {
                queryLabels = {"5", "7", "14"};
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
        else if(lables.size() ==10)
        {
            if(lables.find("15") != lables.end())
            {
                queryLabels = {"0", "4", "10", "14", "5"}; 
            }
            else
            {
                queryLabels = {"0", "2", "4", "6", "8"};
            }
        }else if(lables.size()==12){
            queryLabels = {"0", "2", "4", "6", "8", "10"};
        }else if(lables.size()==11){
            queryLabels = {"2", "3", "9","4"};
        }

        for (int i = 0; i < 5; i++)
        {
            // std::vector<int> validNumbers = hierarchy.vaildCandicate(i, Kmax);

            int k = 0.1 * (i * 2 + 1) * kmax;

            if (k == 1)
                k++;
            else if (k == kmax)
                k--;

            auto query1 = std::chrono::high_resolution_clock::now();

            std::unordered_set<int> solution4 = queryWithLcs(graph1, k, optLcs, queryLabels, core, kmax);

            auto query2 = std::chrono::high_resolution_clock::now();

            // Graph temp(graph1);

            // std::unordered_set<int> solution5 = CommunitySearch::onlineWithEdge(temp, k, queryLabels);

            std::unordered_set<int> solution5 = queryWithMuptIndex(graph1, k, basicLcs, queryLabels);

            auto query3 = std::chrono::high_resolution_clock::now();

            // Graph tempgraph(graph1);

            // std::unordered_set<int> solution6 = CommunitySearch::onlineWithCore(tempgraph, k, queryLabels, core);

            // auto query4 = std::chrono::high_resolution_clock::now();

            if(solution4.size()==solution5.size()){
                indexCost += std::chrono::duration_cast<std::chrono::milliseconds>(query2 - query1).count();
                onlineCost += std::chrono::duration_cast<std::chrono::milliseconds>(query3 - query2).count();
            }
            std::cout << "k为 " << k << " 时,基于索引的查询时间为" << std::chrono::duration_cast<std::chrono::nanoseconds>(query2 - query1).count() << std::endl;
            std::cout << "k为 " << k << " 时,对照算法的查询时间为" << std::chrono::duration_cast<std::chrono::nanoseconds>(query3 - query2).count() << std::endl;
            // std::cout << "k为 " << k << " 时,Baseline的查询时间为" << std::chrono::duration_cast<std::chrono::milliseconds>(query4 - query3).count() << std::endl;
            std::cout << "k为 " << k << " 时,基于索引查询结果数量为" << solution4.size() << std::endl;
            std::cout << "k为 " << k << " 时,对照算法查询结果数量为" << solution5.size() << std::endl;
            // std::cout << "k为 " << k << " 时,Baseline查询结果数量为" << solution6.size() << std::endl;
            std::cout << std::endl;
        }
        std::cout << "基于索引的平均查询时间为" << indexCost / 5 <<" ns"<< std::endl;
        std::cout << "对照算法的平均查询时间为" << onlineCost / 5 <<" ns"<< std::endl;
        std::cout << "平均加速比为" << (onlineCost / 5) / (indexCost / 5) << std::endl;

        return 0;
    }
}