#include "TreeIndex.h"

TreeIndex::TreeIndex(Graph &graph, std::unordered_set<std::string> &labels)
{
    // k-shell分解
    computeCoreIndex(graph, labels);
    // 打印核索引
    // printCoreIndex();
    // printCoreMinimumDegree();

    // 识别并存储每一层的连通分量,这些连通分量属于同一层,相互之间不连通,分量内部连通,且coreIndex/coreMinimumDegree相同
    identifyAndStoreComponents(graph);

    // 打印连通分量数据，以验证正确性
    // printComponents();

    // 构建父子连通分量关系:一个子可能有多个不同的父,不同的父可能有相同的子
    buildParentChildRelationships(graph);

}
void TreeIndex::computeCoreIndex(Graph &graph, std::unordered_set<std::string> &labels)
{
    coreIndex = CoreGroup::coreGroupsAlgorithm(graph);
    // 计算每个核心的最小度数并存储在 coreMinimumDegree 中
    std::set<int> coreIndexSet;
    for (auto &pair : coreIndex)
    {
        coreIndexSet.insert(pair.second);
    }
    int node = 0;
    while (!coreIndexSet.empty())
    {
        coreMinimumDegree[node] = *coreIndexSet.begin();
        for (auto &pair : coreIndex)
        {
            if (pair.second == *coreIndexSet.begin())
            {
                pair.second = node;
            }
        }
        coreIndexSet.erase(coreIndexSet.begin());
        node++;
    }
    //std::cout << "coreMinimumDegree size is :" << coreMinimumDegree.size() << std::endl;

    // 打印最大k
    for (const auto &pair : coreMinimumDegree)
    {
        //std::cout << "核心 " << pair.first << " 的最小度数: " << pair.second << std::endl;
    }
}

void TreeIndex::identifyAndStoreComponents(Graph &graph)
{
    std::unordered_map<int, std::unordered_set<int>> nodesInShell;

    for (const auto &pair : coreIndex)
    {

        nodesInShell[pair.second].insert(pair.first);
    }

    // 根据每一层k-shell的节点，找到每一层的连通分量
    for (const auto &shell : nodesInShell)
    {
        int shellIndex = shell.first;
        const auto &nodes = shell.second;

        // BFS
        std::unordered_set<int> visited;
        for (int node : nodes)
        {

            if (visited.count(node) == 0)
            {
                std::queue<int> queue;
                std::unordered_set<int> component;
                queue.push(node);
                visited.insert(node);

                while (!queue.empty())
                {
                    int currentNode = queue.front();
                    queue.pop();
                    component.insert(currentNode);

                    for (int neighbor : graph.getNeighbors(currentNode))
                    {
                        // 只找未访问过的同一层邻居
                        if (nodes.count(neighbor) > 0 && visited.count(neighbor) == 0)
                        {
                            queue.push(neighbor);
                            visited.insert(neighbor);
                        }
                    }
                }

                int componentId = nextComponentId++;
                layerToComponentToNodes[shellIndex][componentId] = component;
                ComponentToNodes[componentId] = component;
                for (int nodeId : component)
                {
                    nodeToComponentId[nodeId] = componentId;
                }
            }
        }
    }

    //std::cout << "ComponentToNodes size is :" << ComponentToNodes.size() << std::endl;
}

void TreeIndex::printComponents() const
{
    std::cout << "Layer to Component to Nodes Mapping:" << std::endl;
    for (const auto &layer : layerToComponentToNodes)
    {
        std::cout << "Shell Index: " << layer.first << " the k is : " << coreMinimumDegree.at(layer.first) << std::endl;
        for (const auto &comp : layer.second)
        {
            std::cout << "  Component ID: " << comp.first << ", Nodes: ";
            for (int node : comp.second)
            {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
}

void TreeIndex::buildParentChildRelationships(Graph &graph)
{
    // 父子索引的处理
    //std::cout << "Total layers to process: " << layerToComponentToNodes.size() << std::endl;
    int size = layerToComponentToNodes.size();
    std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> lables = graph.getEdgeLables();
    for (int currentLevel = 0; currentLevel < size; ++currentLevel)
    {
        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        // std::cout << "Processing layer " << currentLevel << " with " << currentLayerComponents.size() << " components." << std::endl;

        for (auto &component : currentLayerComponents)
        {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent)
            {
                for (int neighbor : graph.getNeighbors(node))
                {
                    if (coreMinimumDegree[coreIndex[node]] >= coreMinimumDegree[coreIndex[neighbor]])
                    {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId)
                    {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentParent[currentLevel][currentComponentId] = neighborComponentId;
                        ComponentParent[currentComponentId].insert(neighborComponentId);
                        // std::unordered_set<std::string> label = lables[node][neighbor];
                        // parentIndex[node][neighbor].insert(label.begin(),label.end());
                    }
                }
            }
            // 顶层分量的父分量是-1
            if (connectedComponentParent[currentLevel].find(currentComponentId) == connectedComponentParent[currentLevel].end() && ComponentParent[currentComponentId].empty())
            {
                connectedComponentParent[currentLevel][currentComponentId] = -1;
                ComponentParent[currentComponentId].insert(-1);
            }
        }
    }
    for (int currentLevel = size - 1; currentLevel >= 0; currentLevel--)
    {

        auto &currentLayerComponents = layerToComponentToNodes[currentLevel];

        // std::cout << "Processing layer " << currentLevel << " with " << currentLayerComponents.size() << " components." << std::endl;

        for (auto &component : currentLayerComponents)
        {
            int currentComponentId = component.first;
            auto &nodesInCurrentComponent = component.second;

            for (int node : nodesInCurrentComponent)
            {
                // 再打印父索引加入了几个分量 通过这样来判断是否有父子关系遗漏
                for (int neighbor : graph.getNeighbors(node))
                {
                    if (coreMinimumDegree[coreIndex[node]] <= coreMinimumDegree[coreIndex[neighbor]])
                    {
                        continue;
                    }
                    if (nodeToComponentId.find(neighbor) != nodeToComponentId.end() && nodeToComponentId[neighbor] != currentComponentId)
                    {
                        int neighborComponentId = nodeToComponentId[neighbor];
                        connectedComponentChildren[currentLevel][currentComponentId].insert(neighborComponentId);
                        ComponentChildren[currentComponentId].insert(neighborComponentId);
                        // std::unordered_set<std::string> label = lables[node][neighbor];
                        // childIndex[node][neighbor].insert(label.begin(),label.end());
                    }
                }
            }
            // 底层分量的子分量是-1
            if (connectedComponentChildren[currentLevel].find(currentComponentId) == connectedComponentChildren[currentLevel].end())
            {
                connectedComponentChildren[currentLevel][currentComponentId].insert(-1);
                ComponentChildren[currentComponentId].insert(-1);
            }
        }
    }

    //printConnectedComponentChildren();
}

std::unordered_set<int> TreeIndex::findTopComponents(std::vector<int> &queryNodes, int k)
{
    std::unordered_set<int> topComponents;
    std::unordered_set<int> queryComponents;
    std::queue<int> bfsQueue;
    std::unordered_set<int> visited;

    // 初始化：将查询节点的组件ID添加到队列和已访问集合中
    for (int node : queryNodes)
    {
        int componentId = nodeToComponentId[node];
        //std::cout<<"the compId of node is "<<componentId<<std::endl;
        //if (ComponentParent[componentId].find(-1) != ComponentParent[componentId].end())
        {
            // 只有顶层分量（父分量为-1）才加入队列
            queryComponents.insert(componentId);
            bfsQueue.push(componentId);
            visited.insert(componentId);
        }
    }

    while (!bfsQueue.empty())
    {
        int currentComponentId = bfsQueue.front();
        bfsQueue.pop();

        // 探测所有父分量
        for (int parentComponentId : ComponentParent[currentComponentId])
        {
            //std::cout<<"the parent compId is "<<parentComponentId<<std::endl;
            if (parentComponentId == -1)
            {
                // 如果父分量为-1，保存当前分量
                topComponents.insert(currentComponentId);
            }
            if (visited.find(parentComponentId) == visited.end())
            {
                bfsQueue.push(parentComponentId);
                visited.insert(parentComponentId);
            }
        }

        // 探测所有子分量
        for (int childComponentId : ComponentChildren[currentComponentId])
        {
            if (visited.find(childComponentId) == visited.end())
            {
                // 确保所有子分量的节点核心等级必须大于k
                bool eligible = true;
                for (int childNode : ComponentToNodes[childComponentId])
                {
                    if (coreMinimumDegree[coreIndex[childNode]] < k)
                    {
                        eligible = false;
                        break;
                    }
                }
                if (eligible)
                {
                    bfsQueue.push(childComponentId);
                    visited.insert(childComponentId);
                }
            }
        }
    }

    return topComponents;
}

// std::unordered_set<int> TreeIndex::findTopLevelComponents(std::vector<int> &queryNodes, int k)
// {
//     std::unordered_set<int> topComponents;

//     std::unordered_set<int> queryComponents;
//     for (int node : queryNodes)
//     {
//         queryComponents.insert(nodeToComponentId[node]);
//     }

//     // 遍历所有父索引为-1的分量
//     for (const auto &pair : ComponentParent)
//     {
//         int componentId = pair.first;
//         bool isTopLevel = true;

//         // 查看当前分量的父分量是否为-1
//         if (pair.second.find(-1) == pair.second.end())
//         {
//             continue; // 如果父分量不为-1，跳过
//         }

//         // 检查当前分量的所有节点的核心指标是否大于等于 k

//         int node = *ComponentToNodes[componentId].begin();

//         if (coreMinimumDegree[coreIndex[node]] < k)
//         {
//             isTopLevel = false;
//         }
//         else
//         {
//             isTopLevel = true;
//         }

//         // BFS: 查看其子孙的分量是否包括查询节点组包含的分量组
//         std::queue<int> bfsQueue;
//         std::unordered_set<int> visited;

//         // 将子索引加入 BFS 队列
//         for (int childComponentId : ComponentChildren[componentId])
//         {
//             if (ComponentChildren[componentId].find(-1) == ComponentChildren[componentId].end())
//             {
//                 bfsQueue.push(childComponentId);
//             }
//         }

//         while (!bfsQueue.empty())
//         {
//             int currentComponentId = bfsQueue.front();
//             bfsQueue.pop();
//             visited.insert(currentComponentId);

//             // 如果当前分量包含在查询节点组中
//             if (queryComponents.find(currentComponentId) != queryComponents.end())
//             {
//                 isTopLevel = true;
//                 break;
//             }

//             // 将当前分量的子索引加入 BFS 队列
//             for (int childComponentId : ComponentChildren[currentComponentId])
//             {
//                 if (visited.find(childComponentId) == visited.end())
//                 {
//                     bfsQueue.push(childComponentId);
//                 }
//             }
//         }

//         // 如果是顶层分量，则加入结果集合
//         if (isTopLevel)
//         {
//             topComponents.insert(componentId);
//         }
//     }

//     return topComponents;
// }

std::unordered_set<int> TreeIndex::findKCoreSubgraph(std::vector<int> &queryNodes, int k)
{

    int minCoreIndex = INT_MAX;

    // 确定查询节点组的最小核心索引
    for (int node : queryNodes)
    {
        if (coreIndex.find(node) == coreMinimumDegree.end())
        {
            throw std::runtime_error("Error: Node does not exist in core index.");
        }
        minCoreIndex = std::min(minCoreIndex, coreMinimumDegree[coreIndex[node]]);
    }

    if (minCoreIndex < k)
    {
        throw std::runtime_error("Error: The minimum core index of the query nodes is less than the specified k.");
    }

    // 找到所有父分量为-1的顶层分量,开始向下探索孩子分量,查看是否包括某一个的查询节点,通过这样一个个添加顶层分量的节点,最终得到最大k-core子图

    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_set<int> topComponentIds;

    // 向上追溯到顶层
    topComponentIds = findTopComponents(queryNodes, k);

    auto mid = std::chrono::high_resolution_clock::now();

    auto midduration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start);

    //std::cout << "找到父索引所需时间为 :" << midduration.count() << " 毫秒" << std::endl;

    //std::cout << "the topComponentIds size is :" << topComponentIds.size() << std::endl;

    std::unordered_set<int> resultNodes;

    std::unordered_set<int> res;

    std::queue<int> componentQueue;

    // 存储节点id
    for (int compId : topComponentIds)
    {
        componentQueue.push(compId);
    }

    while (!componentQueue.empty())
    {
        int currentId = componentQueue.front();
        componentQueue.pop();

        // 如果当前节点已经在结果集中，则跳过
        if (resultNodes.find(currentId) != resultNodes.end())
            continue;

        // 将当前节点添加到结果集中
        resultNodes.insert(currentId);

        for (int childCompId : ComponentChildren[currentId])
        {
            // 将子连通分量的节点加入队列
            int childNode;
            for (int node : ComponentToNodes[childCompId])
            {
                childNode = node;
                break;
            }

            if (coreMinimumDegree[coreIndex[childNode]] >= k)
            {
                componentQueue.push(childCompId);
            }
        }
    }

    for (int compId : resultNodes)
    {
        for (int node : ComponentToNodes[compId])
        {
            res.insert(node);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();

    // 计算时间差
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid);

    // 输出程序运行时间
    //std::cout << "找到子索引所需时间为 : " << duration.count() << " 毫秒" << std::endl;

    return res;
}

// std::unordered_set<int> TreeIndex::findMaxKCoreSubgraph(std::vector<int> &queryNodes, int k)
// {
//     int minCoreIndex = INT_MAX;

//     // 确定查询节点组的最小核心索引
//     for (int node : queryNodes)
//     {
//         if (coreIndex.find(node) == coreMinimumDegree.end())
//         {
//             throw std::runtime_error("Error: Node does not exist in core index.");
//         }
//         minCoreIndex = std::min(minCoreIndex, coreMinimumDegree[coreIndex[node]]);
//     }

//     if (minCoreIndex < k)
//     {
//         throw std::runtime_error("Error: The minimum core index of the query nodes is less than the specified k.");
//     }

//     // 向上追溯到顶层
//     // 找到所有父分量为-1的顶层分量,开始向下探索孩子分量,查看是否包括某一个的查询节点,通过这样一个个添加顶层分量的节点,最终得到最大k-core子图
//     // TODO
//     auto start = std::chrono::high_resolution_clock::now();
//     std::unordered_set<int> topComponentIds;

//     topComponentIds = findTopComponents(queryNodes, k);

//     auto mid = std::chrono::high_resolution_clock::now();

//     auto midduration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start);

//     std::cout << "找到父索引所需时间为 :" << midduration.count() << " 毫秒" << std::endl;

//     std::cout << "the topComponentIds size is :" << topComponentIds.size() << std::endl;

//     std::unordered_set<int> resultNodes;
//     std::queue<int> componentQueue;

//     // 存储节点id
//     for (int compId : topComponentIds)
//     {
//         for (int node : ComponentToNodes[compId])
//         {
//             componentQueue.push(node);
//         }
//     }

//     while (!componentQueue.empty())
//     {
//         int currentNode = componentQueue.front();
//         componentQueue.pop();

//         // 如果当前节点已经在结果集中，则跳过
//         if (resultNodes.find(currentNode) != resultNodes.end())
//             continue;

//         // 将当前节点添加到结果集中
//         resultNodes.insert(currentNode);

//         // 获取当前节点的子连通分量
//         if (coreMinimumDegree[coreIndex[currentNode]] > k)
//         {
//             // 获取当前节点所在的连通分量ID和层数
//             int currentCompId = nodeToComponentId[currentNode];

//             for (int childCompId : ComponentChildren[currentCompId])
//             {
//                 // 将子连通分量的节点加入队列
//                 for (int childNode : ComponentToNodes[childCompId])
//                 {
//                     if (coreMinimumDegree[coreIndex[childNode]] < k)
//                     {

//                         break;
//                     }
//                     componentQueue.push(childNode);
//                 }
//             }
//         }
//     }
//     auto end = std::chrono::high_resolution_clock::now();

//     // 计算时间差
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid);

//     // 输出程序运行时间
//     std::cout << "找到子索引所需时间为 : " << duration.count() << " 毫秒" << std::endl;

//     return resultNodes;
// }

TreeIndex::TreeIndex(Graph &graph, std::string datasetName)
{
    // computeCoreIndex(graph);

    printCoreIndex();

    computeCoreComposition(graph);

    // printConnectedComponentNodes();
    // printConnectedComponentChildren();
    // printConnectedComponentParent();

    // computeNodeNeighbors(graph);
}

TreeIndex::TreeIndex(const std::string &path)
{
    // deserializeFromFile(path);
}

void TreeIndex::computeCoreComposition(Graph &graph)
{

    std::unordered_map<int, int> highestGroup;

    std::vector<std::stack<int>> coreGroup(graph.getNumberOfNodes());

    for (int i = 0; i < graph.getNumberOfNodes(); ++i)
    {
        coreGroup[i] = std::stack<int>();
    }

    for (const auto &ci : coreIndex)
    {
        coreGroup[coreGroup.size() - ci.second - 1].push(ci.first);
    }

    int connectedComponentIndex = 0;

    int coreIndex = coreGroup.size() - 1;
    for (std::stack<int> &group : coreGroup)
    {
        if (!group.empty())
        {
            connectedComponentChildren[coreIndex] = std::unordered_map<int, std::unordered_set<int>>();
            connectedComponentNodes[coreIndex] = std::unordered_map<int, std::unordered_set<int>>();
            connectedComponentParent[coreIndex] = std::unordered_map<int, int>();
        }

        while (!group.empty())
        {
            int node = group.top();
            // std::cout<<node<<" endl"<<std::endl;
            group.pop();

            for (int neighbor : graph.getNeighbors(node))
            {
                if (this->coreIndex[node] == this->coreIndex[neighbor])
                {
                    if (!nodeGroup.count(neighbor))
                    {
                        continue;
                    }
                    else if (nodeGroup.find(node) == nodeGroup.end())
                    {
                        // Branch 1
                        connectedComponentNodes[coreIndex][nodeGroup[neighbor]].insert(node);
                        nodeGroup[node] = nodeGroup[neighbor];
                        highestGroup[node] = nodeGroup[neighbor];
                    }
                    else if (nodeGroup[node] != nodeGroup[neighbor])
                    {
                        // Branch 2
                        int oldGroup = nodeGroup[neighbor];
                        connectedComponentNodes[coreIndex][nodeGroup[node]].insert(
                            connectedComponentNodes[coreIndex][oldGroup].begin(),
                            connectedComponentNodes[coreIndex][oldGroup].end());
                        connectedComponentChildren[coreIndex][nodeGroup[node]].insert(
                            connectedComponentChildren[coreIndex][oldGroup].begin(),
                            connectedComponentChildren[coreIndex][oldGroup].end());

                        for (int oldGroupNode : connectedComponentNodes[coreIndex][oldGroup])
                        {
                            nodeGroup[oldGroupNode] = nodeGroup[node];
                        }

                        for (auto &oldHighestNode : highestGroup)
                        {
                            if (oldHighestNode.second == oldGroup)
                            {
                                highestGroup[oldHighestNode.first] = nodeGroup[node];
                            }
                        }

                        for (int newParentGroup : this->connectedComponentChildren[coreIndex][oldGroup])
                        {
                            this->connectedComponentParent[coreIndex + 1][newParentGroup] = nodeGroup[node];
                        }

                        connectedComponentNodes[coreIndex].erase(oldGroup);
                        connectedComponentChildren[coreIndex].erase(oldGroup);
                        connectedComponentParent[coreIndex].erase(oldGroup);
                    }
                }
                else
                {
                    if (nodeGroup.find(neighbor) == nodeGroup.end())
                    {
                        continue;
                    }
                    else if (nodeGroup.find(node) == nodeGroup.end())
                    {
                        if (connectedComponentNodes[coreIndex].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex].end())
                        {
                            // Branch 3
                            connectedComponentNodes[coreIndex][highestGroup[neighbor]].insert(node);
                            nodeGroup[node] = highestGroup[neighbor];
                            highestGroup[node] = highestGroup[neighbor];
                        }
                        else if (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex + 1].end())
                        {
                            // Branch 4
                            connectedComponentNodes[coreIndex][connectedComponentIndex].insert(node);
                            connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                            nodeGroup[node] = connectedComponentIndex;
                            highestGroup[node] = connectedComponentIndex;
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = connectedComponentIndex;

                            int oldHighest = highestGroup[neighbor];
                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == oldHighest)
                                {
                                    highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                }
                            }
                            connectedComponentIndex++;
                        }
                        else
                        {
                            // Branch 5

                            int highestLevel = coreIndex + 2;

                            while (connectedComponentNodes[highestLevel].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                highestLevel++; // 不断向上查找，直到找到一个包含邻居最高组的级别
                            }

                            while (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                if (connectedComponentNodes.find(highestLevel - 1) == connectedComponentNodes.end())
                                {
                                    connectedComponentNodes[highestLevel - 1] = std::unordered_map<int, std::unordered_set<int>>();
                                }
                                connectedComponentNodes[highestLevel - 1][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                                connectedComponentParent[highestLevel][highestGroup[neighbor]] = connectedComponentIndex;

                                // 更新highestGroup映射
                                for (auto &oldHighestNode : highestGroup)
                                {
                                    if (oldHighestNode.second == highestGroup[neighbor])
                                    {
                                        highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                    }
                                }

                                highestLevel--;
                                connectedComponentIndex++;
                            }

                            // 创建最初的连接组
                            connectedComponentNodes[coreIndex][connectedComponentIndex] = std::unordered_set<int>({node});
                            connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>({highestGroup[neighbor]});
                            nodeGroup[node] = connectedComponentIndex;
                            highestGroup[node] = connectedComponentIndex;
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = connectedComponentIndex;

                            // 更新highestGroup映射
                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                }
                            }
                            connectedComponentIndex++;
                        }
                    }
                    else if (nodeGroup[node] != highestGroup[neighbor])
                    {
                        if (connectedComponentNodes[coreIndex].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex].end() && nodeGroup[node] != highestGroup[neighbor])
                        {
                            // Branch 6
                            int oldGroup = highestGroup[neighbor];
                            connectedComponentNodes[coreIndex][nodeGroup[node]].insert(
                                connectedComponentNodes[coreIndex][oldGroup].begin(),
                                connectedComponentNodes[coreIndex][oldGroup].end());
                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(
                                connectedComponentChildren[coreIndex][oldGroup].begin(),
                                connectedComponentChildren[coreIndex][oldGroup].end());

                            for (int oldGroupNode : connectedComponentNodes[coreIndex][oldGroup])
                            {
                                nodeGroup[oldGroupNode] = nodeGroup[node];
                            }

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == oldGroup)
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }

                            for (int newParentGroup : this->connectedComponentChildren[coreIndex][oldGroup])
                            {
                                connectedComponentParent[coreIndex + 1][newParentGroup] = nodeGroup[node];
                            }

                            connectedComponentNodes[coreIndex].erase(oldGroup);
                            connectedComponentChildren[coreIndex].erase(oldGroup);
                            connectedComponentParent[coreIndex].erase(oldGroup);
                        }
                        else if (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) != connectedComponentNodes[coreIndex + 1].end())
                        {
                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(highestGroup[neighbor]);
                            connectedComponentParent[coreIndex + 1][highestGroup[neighbor]] = nodeGroup[node];

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }
                        }
                        else
                        { // branch 8
                            int highestLevel = coreIndex + 2;
                            while (connectedComponentNodes[highestLevel].find(highestGroup[neighbor]) == connectedComponentNodes[highestLevel].end())
                            {
                                highestLevel++;
                            }

                            while (connectedComponentNodes[coreIndex + 1].find(highestGroup[neighbor]) == connectedComponentNodes[coreIndex + 1].end())
                            {
                                connectedComponentNodes[highestLevel - 1][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                                connectedComponentChildren[coreIndex][connectedComponentIndex].insert(highestGroup[neighbor]);
                                connectedComponentParent[highestLevel][highestGroup[neighbor]] = connectedComponentIndex;

                                for (auto &oldHighestNode : highestGroup)
                                {
                                    if (oldHighestNode.second == highestGroup[neighbor])
                                    {
                                        highestGroup[oldHighestNode.first] = connectedComponentIndex;
                                    }
                                }

                                highestLevel--;
                                connectedComponentIndex++;
                            }

                            connectedComponentChildren[coreIndex][nodeGroup[node]].insert(highestGroup[neighbor]);
                            connectedComponentParent[coreIndex + 1][highestGroup[node]] = nodeGroup[node];

                            for (auto &oldHighestNode : highestGroup)
                            {
                                if (oldHighestNode.second == highestGroup[neighbor])
                                {
                                    highestGroup[oldHighestNode.first] = nodeGroup[node];
                                }
                            }
                        }
                    }
                }
            }

            // If it has no neighbors
            if (nodeGroup.find(node) == nodeGroup.end())
            {
                connectedComponentNodes[coreIndex][connectedComponentIndex] = std::unordered_set<int>{node};
                connectedComponentChildren[coreIndex][connectedComponentIndex] = std::unordered_set<int>();
                nodeGroup[node] = connectedComponentIndex;
                highestGroup[node] = connectedComponentIndex;
                connectedComponentIndex++;
            }
        }

        coreIndex--;
    }
}

void TreeIndex::computeCoreCompositionByLcy(Graph &graph)
{
    //
}

void printHighestGroup(const std::unordered_map<int, int> &highestGroup)
{
    for (const auto &pair : highestGroup)
    {
        std::cout << "Node: " << pair.first << ", Group: " << pair.second << std::endl;
    }
}

void TreeIndex::computeNodeNeighbors(Graph &graph)
{
    for (auto &node : graph.getNodes())
    {
        for (int neighbor : graph.getNeighbors(node))
        {
            int shellIndex = coreIndex[neighbor];
            nodeNeighbors[node][shellIndex].insert(neighbor);
        }
    }
}

int TreeIndex::getMinimumCoreIndex(std::vector<int> queryNodes)
{
    std::unordered_set<int> queryIndexes;
    for (int queryNode : queryNodes)
    {
        queryIndexes.insert(coreIndex[queryNode]);
    }
    int k = *std::max_element(queryIndexes.begin(), queryIndexes.end());

    std::unordered_set<int> Q(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> C;
    for (int queryNode : queryNodes)
    {
        if (coreIndex[queryNode] == k)
        {
            C.insert(nodeGroup[queryNode]);
            Q.erase(queryNode);
        }
    }

    while (!Q.empty() || C.size() != 1)
    {
        std::unordered_set<int> cFirst;
        for (int connectedComponentID : C)
        {
            cFirst.insert(connectedComponentParent[k][connectedComponentID]);
        }

        for (int queryNode : queryNodes)
        {
            if (coreIndex[queryNode] == k - 1)
            {
                cFirst.insert(nodeGroup[queryNode]);
                Q.erase(queryNode);
            }
        }

        C = cFirst;
        --k;
    }

    return k;
}

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex)
{
    std::unordered_set<int> neighbors;

    for (auto &[shellIndex, nodes] : nodeNeighbors[node])
    {
        if (shellIndex >= coreIndex)
        {
            neighbors.insert(nodes.begin(), nodes.end());
        }
    }

    return neighbors;
}

std::unordered_set<int> TreeIndex::getNeighbors(int node, int coreIndex, std::unordered_set<int> &subcore)
{
    std::unordered_set<int> neighbors;
    for (auto &[shellIndex, nodes] : nodeNeighbors[node])
    {
        if (shellIndex >= coreIndex)
        {
            for (int x : nodes)
            {
                if (subcore.count(x))
                {
                    neighbors.insert(x);
                }
            }
        }
    }

    return neighbors;
}

int TreeIndex::getCoreMinimumDegree(int coreIndex)
{
    return coreMinimumDegree[coreIndex];
}

// 获取 map 中所有的键
std::unordered_set<int> getKeySet(const std::unordered_map<int, std::unordered_set<int>> &map)
{
    std::unordered_set<int> keys;
    for (const auto &pair : map)
    {
        keys.insert(pair.first);
    }
    return keys;
}

// 获取指定核心中的节点数量
int TreeIndex::getNumberOfNodes(int coreIndex)
{
    int currentCoreIndex = coreIndex;
    std::unordered_set<int> currentConnectedComponents = getKeySet(connectedComponentNodes[currentCoreIndex]);

    int numberOfNodes = 0;
    while (!currentConnectedComponents.empty())
    {
        std::unordered_set<int> childrenConnectedComponents;
        for (int connectedComponent : currentConnectedComponents)
        {
            numberOfNodes += connectedComponentNodes[currentCoreIndex][connectedComponent].size();
            if (connectedComponentChildren[currentCoreIndex].count(connectedComponent))
            {
                for (int child : connectedComponentChildren[currentCoreIndex][connectedComponent])
                {
                    childrenConnectedComponents.insert(child);
                }
            }
        }

        currentConnectedComponents = childrenConnectedComponents;
        ++currentCoreIndex;
    }

    return numberOfNodes;
}

std::unordered_set<int> TreeIndex::getCore(std::vector<int> queryNodes)
{
    std::unordered_set<int> core;
    std::unordered_set<int> queryIndexes;
    for (int queryNode : queryNodes)
    {
        queryIndexes.insert(coreIndex[queryNode]);
    }
    int k = *std::max_element(queryIndexes.begin(), queryIndexes.end());

    std::unordered_set<int> Q(queryNodes.begin(), queryNodes.end());
    std::unordered_set<int> C;
    for (int queryNode : queryNodes)
    {
        if (coreIndex[queryNode] == k)
        {
            int cc = nodeGroup[queryNode];
            C.insert(cc);
            Q.erase(queryNode);
        }
    }

    while (!Q.empty() || C.size() != 1)
    {
        std::unordered_set<int> cFirst;
        for (int connectedComponentID : C)
        {
            int parentID = connectedComponentParent[k][connectedComponentID];
            cFirst.insert(parentID);
        }

        for (int queryNode : queryNodes)
        {
            if (coreIndex[queryNode] == k - 1)
            {
                int cc = nodeGroup[queryNode];
                cFirst.insert(cc);
                Q.erase(queryNode);
            }
        }

        C = cFirst;
        --k;
    }

    while (connectedComponentChildren.count(k))
    {
        std::unordered_set<int> newC;
        for (int cc : C)
        {
            newC.insert(cc);
            if (connectedComponentChildren[k].count(cc))
            {
                for (int ccc : connectedComponentChildren[k][cc])
                {
                    newC.insert(ccc);
                }
            }
            if (connectedComponentNodes[k].count(cc))
            {
                core.insert(connectedComponentNodes[k][cc].begin(), connectedComponentNodes[k][cc].end());
            }
        }
        C = newC;
        ++k;
    }

    for (int queryNode : queryNodes)
    {
        core.insert(queryNode);
    }

    return core;
}

int TreeIndex::getCoreIndex(int node)
{
    int index = coreIndex[node];
    return coreMinimumDegree[index];
}

// Get the component ID for a given node
int TreeIndex::getComponent(int node)
{
    auto it = nodeGroup.find(node);
    if (it != nodeGroup.end())
    {
        return it->second;
    }
    else
    {
        return -1; // Or throw an exception or handle error as appropriate
    }
}

// Get the parent component ID for a given component
int TreeIndex::getParentComponent(int componentId)
{
    for (const auto &level : connectedComponentParent)
    {
        auto it = level.second.find(componentId);
        if (it != level.second.end())
        {
            return it->second;
        }
    }
    return 0; // Or handle the case where the component has no parent
}

std::unordered_set<int> TreeIndex::getConnectedComponentChildren(int componentId)
{
    std::unordered_set<int> childrenComponents;

    // 遍历所有核心层
    for (const auto &layer : connectedComponentChildren)
    {
        auto it = layer.second.find(componentId);
        if (it != layer.second.end())
        {
            // 如果找到了给定连通分量ID，加入其所有子连通分量到结果集
            childrenComponents.insert(it->second.begin(), it->second.end());
        }
    }

    return childrenComponents;
}

// Get all nodes in a specified component
std::unordered_set<int> TreeIndex::getNodesInComponent(int componentId)
{
    for (const auto &level : connectedComponentNodes)
    {
        auto it = level.second.find(componentId);
        if (it != level.second.end())
        {
            return it->second;
        }
    }
    return std::unordered_set<int>(); // Return empty if not found
}

// 打印核心指标
void TreeIndex::printCoreIndex()
{
    std::cout << "Core Index:" << std::endl;
    for (const auto &pair : coreIndex)
    {
        std::cout << "Node: " << pair.first << ", Core: " << pair.second << std::endl;
    }
}

// 打印每个核心的最小度数
void TreeIndex::printCoreMinimumDegree()
{
    std::cout << "Core Minimum Degree:" << std::endl;
    for (const auto &pair : coreMinimumDegree)
    {
        std::cout << "Core: " << pair.first << ", Minimum Degree: " << pair.second << std::endl;
    }
}

// 打印连通分量的父子关系
void TreeIndex::printConnectedComponentParent()
{
    std::cout << "Connected Component Parent:" << std::endl;
    for (const auto &outer_pair : connectedComponentParent)
    {
        int parent = outer_pair.first;
        std::cout << "Core Index: " << parent << std::endl;
        const auto &children = outer_pair.second;
        for (const auto &inner_pair : children)
        {
            int child = inner_pair.first;
            std::cout << "    Component ID : " << child << ", Parent: " << inner_pair.second << std::endl;
        }
    }
}

// 打印连通分量的孩子节点集合
void TreeIndex::printConnectedComponentChildren()
{
    std::cout << "Connected Component Children:" << std::endl;
    for (const auto &outer_pair : connectedComponentChildren)
    {
        int parent = outer_pair.first;
        std::cout << "Core Index: " << parent << std::endl;
        const auto &children = outer_pair.second;

        for (const auto &inner_pair : children)
        {
            int child = inner_pair.first;
            const auto &child_set = inner_pair.second;
            std::cout << "    Component ID: " << child << ", Children: ";
            for (int c : child_set)
            {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
}

// 打印连通分量的节点集合
void TreeIndex::printConnectedComponentNodes()
{
    std::cout << "Connected Component Nodes:" << std::endl;
    for (const auto &outer_pair : connectedComponentNodes)
    {
        int connectedComponent = outer_pair.first;
        std::cout << "Core Index: " << connectedComponent << std::endl;
        const auto &nodes = outer_pair.second;
        for (const auto &inner_pair : nodes)
        {
            if (inner_pair.second.size() == 0)
            {
                continue;
            }
            int coreIndex = inner_pair.first;
            const auto &node_set = inner_pair.second;
            std::cout << " Component ID:: " << coreIndex << ", Nodes: ";
            for (int node : node_set)
            {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
}

// 打印节点分组
void TreeIndex::printNodeGroup()
{
    std::cout << "Node Group:" << std::endl;
    for (const auto &pair : nodeGroup)
    {
        std::cout << "Node: " << pair.first << ", Group: " << pair.second << std::endl;
    }
}

// 打印节点邻居信息
void TreeIndex::printNodeNeighbors()
{
    std::cout << "Node Neighbors:" << std::endl;
    for (const auto &outer_pair : nodeNeighbors)
    {
        int node = outer_pair.first;
        std::cout << "Node: " << node << std::endl;
        const auto &core_neighbors = outer_pair.second;
        for (const auto &inner_pair : core_neighbors)
        {
            int coreIndex = inner_pair.first;
            const auto &neighbor_set = inner_pair.second;
            std::cout << "    Core Index: " << coreIndex << ", Neighbors: ";
            for (int neighbor : neighbor_set)
            {
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }
    }
}
void TreeIndex::printQueryNodesCoreIndex(const std::vector<int> &queryNodes)
{
    std::cout << "Core index for query nodes:" << std::endl;
    for (int node : queryNodes)
    {
        std::cout << "Node " << node << ": Core Index = " << coreMinimumDegree.at(coreIndex.at(node)) << " Component is " << nodeToComponentId.at(node) << std::endl;
    }
}
