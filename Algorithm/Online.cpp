#include "Online.h"
// 删除一个节点,检查一次连通性?
// 到达最小度数后,再删除所有不连通的节点

std::unordered_set<int> CommunitySearch::online(Graph &originalGraph, std::vector<int> &queryNodes, int k)
{
    std::unordered_map<int, std::unordered_set<int>> degreeMap;
    std::unordered_map<int, int> nodeDegrees = originalGraph.getDegrees();
    std::unordered_set<int> community = originalGraph.getNodes();

    // 初始化度数映射
    for (int node : community)
    {
        degreeMap[nodeDegrees[node]].insert(node);
    }

    int flag = false;

    int minDegree = 0;

    int nodeToDelete = -1;

    for (int i = 0; i < k; i++)
    {
        if (degreeMap.find(i) != degreeMap.end() && !degreeMap[i].empty())
        {
            flag = true;
            minDegree = i;
            for (int node : degreeMap[minDegree])
            {
                nodeToDelete = node;
                break;
            }
            break;
        }
    }

    // 删除节点直到最小度数达到或超过k
    while (flag)
    {
        // 删除该节点
        degreeMap[minDegree].erase(nodeToDelete);

        community.erase(nodeToDelete);

        // 更新该节点邻居的度数
        for (int neighbor : originalGraph.getNeighbors(nodeToDelete))
        {
            if (community.find(neighbor) != community.end())
            { // 仅更新仍在社区中的节点
                int oldDegree = nodeDegrees[neighbor];
                degreeMap[oldDegree].erase(neighbor);
                nodeDegrees[neighbor]--;
                int newDegree = nodeDegrees[neighbor];
                degreeMap[newDegree].insert(neighbor);
            }
        }
        flag = false;
        for (int i = 0; i < k; i++)
        {
            if (degreeMap.find(i) != degreeMap.end() && !degreeMap[i].empty())
            {
                flag = true;
                minDegree = i;
                for (int node : degreeMap[minDegree])
                {
                    nodeToDelete = node;
                    break;
                }
                break;
            }
        }
    }

    std::queue<int> q;
    std::unordered_set<int> visited;
    std::unordered_set<int> largestComponent;
    for (int node : queryNodes)
    {
        q.push(node);
        visited.insert(node);
    }

    while (!q.empty())
    {
        int current = q.front();
        q.pop();
        largestComponent.insert(current);
        for (int neighbor : originalGraph.getNeighbors(current))
        {
            if (visited.find(neighbor) == visited.end() && community.find(neighbor) != community.end())
            {
                q.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }

    return largestComponent;
}

std::unordered_set<int> CommunitySearch::onlineWithEdge(Graph &originalGraph, int k, std::unordered_set<std::string> &lables)
{

    originalGraph.filterEdgeLableAndRecover(lables);

    // 此时图内只剩下符合边标签集约束的边

    // 得到度数-节点集合 节点-度数集合 节点集合
    std::vector<std::unordered_set<int>> degreesVec = originalGraph.getOrderedNodes();
    std::unordered_map<int, int> nodeDegrees = originalGraph.getDegrees();
    std::unordered_set<int> community = originalGraph.getNodes();

    int flag = false;

    int minDegree = 0;

    int nodeToDelete = -1;

    // 找到第一个待删除的节点
    for (int i = 0; i <= k; i++)
    {
        if (i == k)
        {
            return {};
        }
        if (!degreesVec[i].empty())
        {
            flag = true;
            minDegree = i;
            for (int node : degreesVec[minDegree])
            {
                nodeToDelete = node;
                break;
            }
            break;
        }
    }

    // 删除节点直到最小度数达到或超过k
    while (flag)
    {
        // std::cout << community.size() << std::endl;
        // std::cout<<"the node to delete is "<<nodeToDelete<<std::endl;
        community.erase(nodeToDelete);
        degreesVec[nodeDegrees[nodeToDelete]].erase(nodeToDelete);
        for (int neighbor : originalGraph.getNeighbors(nodeToDelete))
        {
            if (community.find(neighbor) != community.end())
            { // 仅更新仍在社区中的节点
                degreesVec[nodeDegrees[neighbor]].erase(neighbor);
                nodeDegrees[neighbor]--;
                degreesVec[nodeDegrees[neighbor]].insert(neighbor);
            }
        }
        flag = false;
        for (int i = 0; i < k; i++)
        {
            if (!degreesVec[i].empty())
            {
                flag = true;
                minDegree = i;
                for (int node : degreesVec[minDegree])
                {
                    nodeToDelete = node;
                    break;
                }
                break;
            }
        }
    }

    return community;
}

std::unordered_set<int> CommunitySearch::onlineWithCore(Graph &originalGraph, int k, std::unordered_set<std::string> &lables, std::vector<std::unordered_set<int>> &coreness)
{
    originalGraph.filterEdgeLableAndRecover(lables);
    std::unordered_set<int> community;
    // 此时图内只剩下符合边标签集约束的边
    for (int i = k; i < coreness.size(); i++)
    {
        if (!coreness[i].empty())
        {
            community.insert(coreness[i].begin(), coreness[i].end());
        }
    }

    std::vector<std::unordered_set<int>> degreesVec(originalGraph.getDegree() + 1);
    std::unordered_map<int, int> nodeDegrees;
    std::unordered_map<int, std::unordered_set<int>> &adj = originalGraph.getAdj();
    for (int node : community)
    {
        int degree = 0;
        for (int neighbor : adj[node])
        {
            if (community.find(neighbor) != community.end())
            {
                degree++;
            }
        }

        {
            degreesVec[degree].insert(node);
            nodeDegrees[node] = degree;
        }
    }

    int flag = false;

    int minDegree = 0;

    int nodeToDelete = -1;

    // 找到第一个待删除的节点
    for (int i = 0; i <= k; i++)
    {
        if (i == k)
        {
            return {};
        }
        if (!degreesVec[i].empty())
        {
            flag = true;
            minDegree = i;
            for (int node : degreesVec[minDegree])
            {
                nodeToDelete = node;
                break;
            }
            break;
        }
    }

    // 删除节点直到最小度数达到或超过k
    while (flag)
    {
        // std::cout << community.size() << std::endl;
        // std::cout<<"the node to delete is "<<nodeToDelete<<std::endl;
        community.erase(nodeToDelete);
        degreesVec[nodeDegrees[nodeToDelete]].erase(nodeToDelete);
        for (int neighbor : originalGraph.getNeighbors(nodeToDelete))
        {
            if (community.find(neighbor) != community.end())
            { // 仅更新仍在社区中的节点
                degreesVec[nodeDegrees[neighbor]].erase(neighbor);
                nodeDegrees[neighbor]--;
                degreesVec[nodeDegrees[neighbor]].insert(neighbor);
            }
        }
        flag = false;
        for (int i = 0; i < k; i++)
        {
            if (!degreesVec[i].empty())
            {
                flag = true;
                minDegree = i;
                for (int node : degreesVec[minDegree])
                {
                    nodeToDelete = node;
                    break;
                }
                break;
            }
        }
    }

    return community;
}

/*
首先删除所有不满足lables约束的边，然后按照online的方法删除节点，最后返回最大连通子图
*/
std::unordered_set<int> CommunitySearch::onlineWithEdge(Graph &originalGraph, std::vector<int> &queryNodes, int k, std::unordered_set<std::string> &lables)
{

    originalGraph.filterEdgeLableAndRecover(lables);

    // 此时图内只剩下符合边标签集约束的边

    // 得到度数-节点集合 节点-度数集合 节点集合
    std::vector<std::unordered_set<int>> degreesVec = originalGraph.getOrderedNodes();
    std::unordered_map<int, int> nodeDegrees = originalGraph.getDegrees();
    std::unordered_set<int> community = originalGraph.getNodes();

    int flag = false;

    int minDegree = 0;

    int nodeToDelete = -1;

    // 找到第一个待删除的节点
    for (int i = 0; i < k; i++)
    {
        if (!degreesVec[i].empty())
        {
            flag = true;
            minDegree = i;
            for (int node : degreesVec[minDegree])
            {
                nodeToDelete = node;
                break;
            }
            break;
        }
    }

    // 删除节点直到最小度数达到或超过k
    while (flag)
    {
        // std::cout << community.size() << std::endl;
        // std::cout<<"the node to delete is "<<nodeToDelete<<std::endl;
        community.erase(nodeToDelete);
        degreesVec[nodeDegrees[nodeToDelete]].erase(nodeToDelete);
        for (int neighbor : originalGraph.getNeighbors(nodeToDelete))
        {
            if (community.find(neighbor) != community.end())
            { // 仅更新仍在社区中的节点
                degreesVec[nodeDegrees[neighbor]].erase(neighbor);
                nodeDegrees[neighbor]--;
                degreesVec[nodeDegrees[neighbor]].insert(neighbor);
            }
        }
        flag = false;
        for (int i = 0; i < k; i++)
        {
            if (!degreesVec[i].empty())
            {
                flag = true;
                minDegree = i;
                for (int node : degreesVec[minDegree])
                {
                    nodeToDelete = node;
                    break;
                }
                break;
            }
        }
    }

    if (community.find(*queryNodes.begin()) == community.end())
    {
        return {};
    }

    // 从查询节点开始检索最大连通分量
    std::queue<int> q;
    std::unordered_set<int> visited;
    std::unordered_set<int> largestComponent;
    for (int node : queryNodes)
    {
        q.push(node);
        visited.insert(node);
    }

    while (!q.empty())
    {
        // std::cout << "BFS" << std::endl;
        int current = q.front();
        q.pop();
        largestComponent.insert(current);
        for (int neighbor : originalGraph.getNeighbors(current))
        {
            if (visited.find(neighbor) == visited.end() && community.find(neighbor) != community.end())
            {
                q.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }

    return largestComponent;
}

bool CommunitySearch::isConnected(Graph &graph, const std::vector<int> &queryNodes, const std::unordered_set<int> &subgraphNodes)
{
    // std::cout<<"isConnected?"<<std::endl;
    if (queryNodes.empty() || subgraphNodes.empty())
    {
        return false;
    }

    std::set<int> visited;
    std::queue<int> queue;
    for (int startNode : queryNodes)
    {
        if (subgraphNodes.find(startNode) != subgraphNodes.end())
        {
            queue.push(startNode);
            visited.insert(startNode);
            break;
        }
    }

    if (queue.empty())
    {
        return false; // 如果没有有效的起始节点，则直接返回不连通
    }

    while (!queue.empty())
    {
        // std::cout<<"BFS"<<std::endl;
        int node = queue.front();
        queue.pop();
        for (int neighbor : graph.getNeighbors(node))
        {
            if (subgraphNodes.find(neighbor) != subgraphNodes.end() && visited.insert(neighbor).second)
            {
                queue.push(neighbor);
            }
        }
    }

    // 检查所有查询节点是否都在访问集合中
    for (int node : queryNodes)
    {
        if (visited.find(node) == visited.end())
        {
            return false;
        }
    }

    return true;
}
