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
    std::string opt_flag = argv[2];
    std::string index_flag = argv[3];

    int edge = 2689454;

    Graph graph1(dataset_path,edge*0.2);
    Graph graph2(dataset_path,edge*0.4);
    Graph graph3(dataset_path,edge*0.6);
    Graph graph4(dataset_path,edge*0.8);

    auto lcsStart = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::string> labels = graph1.checkLabelForm();

    // 创建和销毁对象，并测量时间
    auto createAndMeasure = [](const std::unordered_set<std::string>& lbls, Graph& g, int mode) {
        auto start = std::chrono::high_resolution_clock::now();
        LCS* lcs;
        if(mode==0||mode==1) 
            lcs = new LCS(lbls, g, mode);
        else 
            lcs = new LCS(lbls, g);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = end - start;
        std::cout << "Time to create LCS: " << time_span.count() << " milliseconds." << std::endl;
        delete lcs;
    };

    //createAndMeasure(labels, graph1, 0); // 对应于 Lcs1
    //createAndMeasure(labels, graph2, 0); // 对应于 Lcs2
    //createAndMeasure(labels, graph3, 0); // 对应于 Lcs3
    //createAndMeasure(labels, graph4, 0); // 对应于 Lcs4

    //createAndMeasure(labels, graph1, 2); // 对应于 Lcs1
    //createAndMeasure(labels, graph2, 2); // 对应于 Lcs2
    createAndMeasure(labels, graph3, 2); // 对应于 Lcs3
    //createAndMeasure(labels, graph4, 2); // 对应于 Lcs4



    return 0;

    // std::cout << "the basic LCS Index built in " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart).count() << " ms" << std::endl;
    // size_t memoryUsage = optLcs1.calculateMemoryUsage();
    // std::cout << "Memory usage of LcsIndex: " << memoryUsage / (1024.0 * 1024.0) << " MB" << std::endl;
    // memoryUsage = optLcs2.calculateMemoryUsage();
    // std::cout << "Memory usage of LcsIndex: " << memoryUsage / (1024.0 * 1024.0) << " MB" << std::endl;
    // memoryUsage = optLcs3.calculateMemoryUsage();
    // std::cout << "Memory usage of LcsIndex: " << memoryUsage / (1024.0 * 1024.0) << " MB" << std::endl;
    // memoryUsage = optLcs4.calculateMemoryUsage();
    // std::cout << "Memory usage of LcsIndex: " << memoryUsage / (1024.0 * 1024.0) << " MB" << std::endl;

}