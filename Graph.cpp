#include "Graph.h"

Graph::Graph(std::string path)
{

    readFromFile(path);
    computeDegrees();
    computeMinimumDegree();
    statistic();
}

Graph::Graph(std::string path,int radio)
{

    readFromFile(path,radio);
    computeDegrees();
    computeMinimumDegree();
    statistic();
}

Graph::Graph(std::string path,int v,int u,std::string label)
{
    readFromFile(path,v,u,label);
    computeDegrees();
    computeMinimumDegree();
    statistic();
}

std::unordered_map<int, std::unordered_set<int>>& Graph::getAdj(){
    return adj;
}

void Graph::statistic()
{
    // 节点数
    std::cout << "节点数为 " << adj.size() << std::endl;
    // 边数
    int edges = 0;
    for (auto &entry : labelsForEdges)
    {
        // for(auto &pair : entry.second){
        //     edges += pair.second.size();
        // }
        edges += entry.second.size();
    }
    std::cout << "边数为 " << edges / 2 << std::endl;

    m = edges / 2;
    // 单个标签最多边数
    int max = 0;
    int avg = 0;
    for (auto &entry : labelsForEdges)
    {
        if (entry.second.size() > max)
        {
            max = entry.second.size();
        }
        avg += entry.second.size();
        std::cout << "标签 " << entry.first << " 的边数量为 " << entry.second.size() / 2 << std::endl;
    }
    std::cout << "单个标签最多边数为 " << max / 2 << std::endl;
    // 单个标签平均边数
    avg /= labelsForEdges.size();
    std::cout << "标签数为 " << labelsForEdges.size() << std::endl;
    std::cout << "单个标签平均边数为 " << avg / 2 << std::endl;
    // 去除多边的边数
    int sum = 0;
    for (auto &entry : adj)
    {
        sum += entry.second.size();
    }
    std::cout << "去除多边的边数为 " << sum / 2 << std::endl;
    // 最大度数
    std::cout << "最大度数为 " << Dmax << std::endl;
}

int Graph::getEdgeNum()
{
    return m;
}

// 复制构造函数的实现
Graph::Graph(Graph &graph)
{
    // 拷贝邻接表
    adj = graph.adj; // unordered_map 自动进行深拷贝

    // 拷贝节点的度
    degrees = graph.degrees; // unordered_map 的深拷贝

    // 拷贝节点标签
    nodeLables = graph.nodeLables; // 深拷贝unordered_map

    // 拷贝边标签
    edgeLables = graph.edgeLables; // 深拷贝内部unordered_map和集合

    // 拷贝有相同度的节点的集合的向量
    orderedNodes = graph.orderedNodes; // 深拷贝向量和集合

    // 拷贝最小度
    minimumDegree = graph.minimumDegree;

}

Graph::~Graph() {}

void Graph::readFromFile(std::string fileName,int radio)
{
    // std::cout << "Begin to Read Data File" << std::endl;
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    size_t total_memory = 0;
    n = 0;
    std::string line;
    // 测试数据集的边是否适用无向图
    while (std::getline(file, line))
    {
        sum++;
        if(radio<=sum)break;
        int from, to;
        char buffer[256]; // Assuming the label will not be longer than 255 characters
        if (sscanf(line.c_str(), "%d %d %255s", &from, &to, buffer) == 3)
        {
            if (0 <= from && 0 <= to && from != to)
            {
                ;
            }
            else
            {
                continue;
            }
            std::string lable(buffer); // Convert char array to std::string
            size_t pos = lable.find("to");
            if (pos != std::string::npos)
            {

                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);
                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge);
                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge);
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
            else
            {
                // Handle the case where "to" is not found in the label
                // std::cout << " Label does not contain 'to'" << std::endl;
                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);

                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge); // 插入边 from -> to

                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge); // 插入边 to -> from
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
        }
        else
        {
            // Handle the case where input format does not match
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    // std::cout << sum << std::endl;
    // std::cout << Double << std::endl;
    // std::cout << "the size of nodes is " << nodes.size() << std::endl;
    // std::cout << "the size of edges is " << sum << std::endl;
    file.close();
}

void Graph::readFromFile(std::string fileName)
{
    // std::cout << "Begin to Read Data File" << std::endl;
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    n = 0;
    std::string line;
    // 测试数据集的边是否适用无向图
    while (std::getline(file, line))
    {
        sum++;
        int from, to;
        char buffer[256]; // Assuming the label will not be longer than 255 characters
        if (sscanf(line.c_str(), "%d %d %255s", &from, &to, buffer) == 3)
        {
            if (0 <= from && 0 <= to && from != to)
            {
                ;
            }
            else
            {
                continue;
            }
            std::string lable(buffer); // Convert char array to std::string
            size_t pos = lable.find("to");
            if (pos != std::string::npos)
            {

                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);
                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge);
                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge);
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
            else
            {
                // Handle the case where "to" is not found in the label
                // std::cout << " Label does not contain 'to'" << std::endl;
                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);

                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge); // 插入边 from -> to

                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge); // 插入边 to -> from
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
        }
        else
        {
            // Handle the case where input format does not match
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    // std::cout << sum << std::endl;
    // std::cout << Double << std::endl;
    // std::cout << "the size of nodes is " << nodes.size() << std::endl;
    // std::cout << "the size of edges is " << sum << std::endl;
    file.close();
}

void Graph::readFromFile(std::string fileName,int u,int v,std::string label)
{
    // std::cout << "Begin to Read Data File" << std::endl;
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error: File not found." << std::endl;
        exit(1);
    }
    std::unordered_set<int> nodes;
    int sum = 0;
    n = 0;
    std::string line;
    //添加增量边
    {
        addNode(u);
        addNode(v);

        addEdge(u,v);
        addEdgeLable(u, v, label);

        addEdge(v,u);
        addEdgeLable(v, u, label);

        nodes.insert(u);
        nodes.insert(v);
        std::pair<int, int> edge(u, v);
        labelsForEdges[label].push_back(edge);
        std::pair<int, int> reverseEdge(v, u);
        labelsForEdges[label].push_back(reverseEdge);
    }
    // 测试数据集的边是否适用无向图
    while (std::getline(file, line))
    {
        sum++;
        int from, to;
        char buffer[256]; // Assuming the label will not be longer than 255 characters
        if (sscanf(line.c_str(), "%d %d %255s", &from, &to, buffer) == 3)
        {
            if (0 <= from && 0 <= to && from != to)
            {
                ;
            }
            else
            {
                continue;
            }
            std::string lable(buffer); // Convert char array to std::string
            size_t pos = lable.find("to");
            if (pos != std::string::npos)
            {

                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);
                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge);
                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge);
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
            else
            {
                // Handle the case where "to" is not found in the label
                // std::cout << " Label does not contain 'to'" << std::endl;
                addNode(from);

                addNode(to);

                addEdge(from, to);
                addEdgeLable(from, to, lable);

                addEdge(to, from);
                addEdgeLable(to, from, lable);

                nodes.insert(from);
                nodes.insert(to);

                std::pair<int, int> edge(from, to);
                labelsForEdges[lable].push_back(edge); // 插入边 from -> to

                std::pair<int, int> reverseEdge(to, from);
                labelsForEdges[lable].push_back(reverseEdge); // 插入边 to -> from
                if (from > n)
                    n = from;
                if (to > n)
                    n = to;
                if (from > to)
                {
                    std::pair<int, int> reverseEdge(to, from);
                    labelsForDirectedEdges[lable].push_back(reverseEdge);
                }
                else
                {
                    std::pair<int, int> edge(from, to);
                    labelsForDirectedEdges[lable].push_back(edge);
                }
            }
        }
        else
        {
            // Handle the case where input format does not match
            std::cerr << "Error: Incorrect line format" << std::endl;
        }
    }
    // std::cout << sum << std::endl;
    // std::cout << Double << std::endl;
    // std::cout << "the size of nodes is " << nodes.size() << std::endl;
    // std::cout << "the size of edges is " << sum << std::endl;
    file.close();
}

bool Graph::addEdge(int from, int to)
{
    if (from != to && adj[from].find(to) == adj[from].end())
    {
        adj[from].insert(to);
        return true;
    }else return false;
}

std::unordered_map<std::string, std::vector<std::pair<int, int>>> &Graph::getLabelsForDirectedEdges()
{
    return labelsForDirectedEdges;
}

void Graph::addNode(int node)
{
    if (adj.find(node) == adj.end())
    {
        adj[node] = std::unordered_set<int>();
    }
}

void Graph::printGraphInfo(Graph &graph)
{
    long long edgeSize = 0;
    for (auto &entry : edgeLables)
    {
        for (auto &ptr : entry.second)
        {
            edgeSize += ptr.second.size();
        }
    }
    edgeSize /= 2;
    std::cout << "开始输出图" << std::endl;
    // for (auto ptr : edges)
    // {
    //     std::cout << "Edge: " << ptr.first << " -> " << ptr.second << std::endl;
    // }
    // for (auto node : graph.getNodes())
    // {
    //     std::cout << "Node : " << node << std::endl;
    // }
    std::cout << "图中边的数量: " << edgeSize << std::endl;
    std::cout << "图中节点的数量: " << graph.getNodes().size() << std::endl;
}

// TODO
std::vector<std::pair<int, int>> Graph::getEdges()
{
    std::vector<std::pair<int, int>> edges;
    std::unordered_set<int> printedNodes;

    for (const auto &entry : adj)
    {
        int from = entry.first;
        for (int neighbor : entry.second)
        {
            if (printedNodes.find(neighbor) == printedNodes.end())
            {
                edges.push_back(std::make_pair(from, neighbor));
            }
        }
        printedNodes.insert(entry.first);
    }

    return edges;
}

void Graph::filteEdges(std::string label)
{
    std::vector<std::pair<int, int>> edges;
    std::unordered_set<int> printedNodes;

    for (const auto &entry : adj)
    {
        int from = entry.first;
        for (int neighbor : entry.second)
        {
            if (printedNodes.find(neighbor) == printedNodes.end())
            {
                edges.push_back(std::make_pair(from, neighbor));
            }
        }
        printedNodes.insert(entry.first);
    }
}

bool Graph::isExistEdge(int v1, int v2)
{
    if (adj[v1].find(v2) != adj[v1].end() && adj[v2].find(v1) != adj[v2].end())
        return true;
    else
        return false;
}

// Add Labels
std::unordered_map<int, int> Graph::computeDegrees()
{
    // std::cout << "Begin to Compute Degrees" << std::endl;
    degrees.clear();
    orderedNodes.clear();
    orderedNodes.resize(adj.size());

    Dmax = 0;

    for (auto &entry : adj)
    {
        degrees[entry.first] = entry.second.size();
        orderedNodes[entry.second.size()].insert(entry.first);
        if (entry.second.size() > Dmax)
            Dmax = entry.second.size();
    }

    return degrees;
}

int Graph::computeMinimumDegree()
{
    // std::cout << "Begin to Compute min Degeree" << std::endl;
    for (size_t i = 0; i < orderedNodes.size(); ++i)
    {
        if (!orderedNodes[i].empty())
        {
            minimumDegree = i;
            return i;
        }
    }
    return minimumDegree;
}

bool Graph::doesNodeExist(int node)
{
    return adj.find(node) != adj.end();
}

std::unordered_map<int, int> Graph::getDegrees()
{
    return degrees;
}

int Graph::getDegree(int node)
{
    return adj[node].size();
}

int Graph::getMinimumDegree()
{
    return minimumDegree;
}

std::unordered_set<int> Graph::getNeighbors(int node)
{
    // 获取并返回排序后的邻居节点
    //std::vector<int> sortedNeighbors = sortNeighbors(node);
    return adj[node];
}

std::vector<int> Graph::sortNeighbors(int node)
{
    // 获取节点的邻居节点
    std::unordered_set<int> &neighbors = adj[node];

    // 将邻居节点存储到向量中
    std::vector<int> sortedNeighbors(neighbors.begin(), neighbors.end());

    // 对向量中的元素进行排序
    std::sort(sortedNeighbors.begin(), sortedNeighbors.end());

    return sortedNeighbors;
}

int Graph::getNumberOfNodes()
{
    return adj.size();
}

int Graph::getNumberOfEdges()
{
    int numberOfEdges = 0;
    for (auto &entry : adj)
    {
        int sum = 0;
        for (auto &ptr : entry.second)
        {
            sum += edgeLables[entry.first][ptr].size();
        }
        numberOfEdges += sum;
    }
    return numberOfEdges / 2;
}

// 寻找图中的最大连通分量
std::unordered_set<int> Graph::findLargestConnectedComponent()
{
    std::unordered_set<int> largestComponent;
    std::unordered_set<int> visited;
    std::queue<int> queue;

    // 遍历图中的每个节点，找到最大的连通分量
    for (const auto &entry : adj)
    {
        if (visited.find(entry.first) == visited.end())
        {
            std::unordered_set<int> component;
            queue.push(entry.first);

            while (!queue.empty())
            {
                int current = queue.front();
                queue.pop();

                if (visited.find(current) == visited.end())
                {
                    visited.insert(current);
                    component.insert(current);

                    // 将所有未访问的邻居加入队列
                    for (int neighbor : adj[current])
                    {
                        if (visited.find(neighbor) == visited.end())
                        {
                            queue.push(neighbor);
                        }
                    }
                }
            }

            // 检查这个组件是否是最大的
            if (component.size() > largestComponent.size())
            {
                largestComponent = std::move(component);
            }
        }
    }

    return largestComponent;
}

int Graph::getNumberOfEdgesOfConnectedComponent(std::unordered_set<int> connectedComponent)
{
    int numberOfEdges = 0;
    for (auto &entry : adj)
    {
        if (connectedComponent.find(entry.first) != connectedComponent.end())
        {
            int sum = 0;
            for (auto &ptr : entry.second)
            {
                sum += edgeLables[entry.first][ptr].size();
            }
            numberOfEdges += sum;
        }
    }
    return numberOfEdges / 2;
}

std::unordered_set<int> Graph::getNodes()
{
    std::unordered_set<int> nodes;
    for (auto &entry : adj)
    {
        nodes.insert(entry.first);
    }
    return nodes;
}

void Graph::removeNode(int node)
{
    // 如果节点不存在于图中，则直接返回
    if (adj.find(node) == adj.end())
    {
        return;
    }

    // 删除与节点相关的所有边
    for (int neighbor : adj[node])
    {
        adj[neighbor].erase(node); // 从邻居节点的邻接列表中删除该节点
        orderedNodes[degrees[neighbor]].erase(neighbor);
        degrees[neighbor]--;
        orderedNodes[degrees[neighbor]].insert(neighbor);
        //edgeLables[neighbor].erase(node);
    }

    // 删除节点及其相关边
    adj.erase(node);
    orderedNodes[degrees[node]].erase(node);
    degrees.erase(node);
}

// 删单向边
void Graph::removeEdge(int from, int to, std::string label)
{
    if (adj[from].find(to) != adj[from].end())
    {
        // 删除标签，检查是否需要完全移除边
        auto &labels = edgeLables[from][to];
        labels.erase(label);

        if (labels.empty())
        {
            // 减少to的度数并更新有序节点集
            orderedNodes[degrees[from]].erase(from);
            degrees[from]--;
            orderedNodes[degrees[from]].insert(from);

            adj[from].erase(to);
            edgeLables[from].erase(to);
        }
        else
        {
            edgeLables[from][to].erase(label);
        }
    }
}

std::vector<std::unordered_set<int>> Graph::getOrderedNodes()
{
    return this->orderedNodes;
}

int Graph::getDegree(){
    return orderedNodes.size();
}

void Graph::serializeToFile(std::string datasetName)
{
    std::ofstream file(datasetName + "Dataset.bin", std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error: Failed to open file." << std::endl;
        exit(1);
    }

    for (auto &entry : adj)
    {
        file.write(reinterpret_cast<const char *>(&entry.first), sizeof(int));
        int size = entry.second.size();
        file.write(reinterpret_cast<const char *>(&size), sizeof(int));
        for (int neighbor : entry.second)
        {
            file.write(reinterpret_cast<const char *>(&neighbor), sizeof(int));
        }
    }

    file.close();
}

void Graph::printDegrees()
{
    std::cout << "Node degrees:" << std::endl;
    for (const auto &pair : degrees)
    {
        std::cout << "Node " << pair.first << ": Degree " << pair.second << std::endl;
    }
}

void Graph::printOrderedNodes()
{
    std::cout << "Ordered nodes:" << std::endl;
    for (size_t i = 0; i < orderedNodes.size(); ++i)
    {
        std::cout << "Nodes with degree " << i << ": ";
        for (int node : orderedNodes[i])
        {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}

// 最好不用
void Graph::keepLargestComponent()
{
    std::vector<std::unordered_set<int>> components = computeConnectedComponents();
    if (components.size() > 1)
    {
        int maxSize = 0;
        std::unordered_set<int> *largestComponent = nullptr;
        for (auto &component : components)
        {
            if (component.size() > maxSize)
            {
                maxSize = component.size();
                if (largestComponent)
                {
                    delete largestComponent; // 删除之前的最大组件
                }
                largestComponent = new std::unordered_set<int>(component); // 分配新的内存
            }
        }
        if (largestComponent)
        { // 确保 largestComponent 不为空
            // 从图中删除非最大连通分量中的节点
            for (int node : getNodes())
            {
                if (largestComponent->find(node) == largestComponent->end())
                {
                    removeNode(node);
                }
            }
            delete largestComponent; // 删除最大组件的内存
        }
    }
}

// 最好不用,因为有dfs可能爆栈
std::vector<std::unordered_set<int>> Graph::computeConnectedComponents()
{
    std::vector<std::unordered_set<int>> components;
    std::unordered_set<int> visited;
    for (int node : getNodes())
    {
        if (visited.find(node) == visited.end())
        {
            std::unordered_set<int> *component = new std::unordered_set<int>();
            dfs(node, visited, *component);
            components.push_back(*component);
            delete component; // 在不再需要时释放内存
        }
    }
    return components;
}

bool Graph::isGraphConnected()
{
    std::unordered_set<int> visited;
    dfsStack(*getNodes().begin(), visited);
    std::cout << "the size of visited is " << visited.size() << std::endl;
    return visited.size() == getNumberOfNodes();
}

void Graph::dfsStack(int startNode, std::unordered_set<int> &visited)
{
    std::stack<int> nodeStack;
    nodeStack.push(startNode);

    while (!nodeStack.empty())
    {
        int currentNode = nodeStack.top();
        nodeStack.pop();

        if (visited.find(currentNode) == visited.end())
        {
            visited.insert(currentNode);

            for (int neighbor : getNeighbors(currentNode))
            {
                if (visited.find(neighbor) == visited.end())
                {
                    nodeStack.push(neighbor);
                }
            }
        }
    }
}

// 最好不用
void Graph::dfs(int node, std::unordered_set<int> &visited, std::unordered_set<int> &component)
{
    visited.insert(node);
    component.insert(node);
    for (int neighbor : getNeighbors(node))
    {
        if (visited.find(neighbor) == visited.end())
        {
            dfs(neighbor, visited, component);
        }
    }
}

std::unordered_set<std::string> Graph::getNodesLable(int id)
{
    return nodeLables[id];
}

std::unordered_set<std::string> Graph::getEdgeLable(int from, int to)
{
    return edgeLables[from][to];
}

void Graph::addEdgeLable(int from, int to, std::string lable)
{
    edgeLables[from][to].insert(lable);
}

void Graph::addNodeLable(int id, std::string lable)
{
    nodeLables[id].insert(lable);
}

// 删除所有不符合边标签集约束的边
void Graph::filterEdgeLable(std::unordered_set<std::string> lables)
{
    std::vector<std::tuple<int, int, std::string>> tuple_list;

    for (int node : getNodes())
    {
        for (int neighbor : getNeighbors(node))
        {
            // std::cout<<"arrive"<<std::endl;
            std::unordered_set<std::string> edgeLables = getEdgeLable(node, neighbor);
            bool flag = false;
            for (std::string lable : edgeLables)
            {
                if (lables.find(lable) == lables.end())
                {
                    // graph.removeEdge(node, neighbor, lable);
                    tuple_list.push_back(std::make_tuple(node, neighbor, lable));
                }
                else
                {
                    flag = true;
                }
            }
        }
    }
    for (auto tuple : tuple_list)
    {
        removeEdge(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple));
        removeEdge(std::get<1>(tuple), std::get<0>(tuple), std::get<2>(tuple));
    }
} // namespace CoreGroup

std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &Graph::getEdgeLables()
{
    return edgeLables;
}
// 删除所有不符合边标签集约束的边
void Graph::filterEdgeLableAndRecover(std::unordered_set<std::string> lables)
{
    std::vector<std::tuple<int, int, std::string>> tuple_list;

    int edges = 0;

    for (auto &entry : edgeLables)
    {
        edges += entry.second.size();
    }
    // std::cout << "before filter , the edge's size is " << edges << std::endl;

    for (int node : getNodes())
    {
        for (int neighbor : getNeighbors(node))
        {
            // std::cout<<"arrive"<<std::endl;
            std::unordered_set<std::string> edgeLables = getEdgeLable(node, neighbor);
            bool flag = false;
            for (std::string lable : edgeLables)
            {
                if (lables.find(lable) == lables.end())
                {
                    // graph.removeEdge(node, neighbor, lable);
                    tuple_list.push_back(std::make_tuple(node, neighbor, lable));
                }
                else
                {
                    flag = true;
                }
            }
        }
    }
    for (auto tuple : tuple_list)
    {
        removeEdge(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple));
        removeEdge(std::get<1>(tuple), std::get<0>(tuple), std::get<2>(tuple));
    }

    edges = 0;

    for (auto &entry : edgeLables)
    {
        edges += entry.second.size();
    }
    // std::cout << "after filter , the edge's size is " << edges << std::endl;
}

std::unordered_map<std::string, std::vector<std::pair<int, int>>> Graph::getLabelsForEdges()
{
    return labelsForEdges;
}

std::vector<std::pair<int, int>> Graph::addEdgeSet(std::vector<std::pair<int, int>> edgeSet, std::string label)
{
    std::vector<std::pair<int, int>> edges;
    for (auto &entry : edgeSet)
    {
        if (addEdge(entry.first, entry.second))
        {
            const auto e = std::make_pair(entry.first, entry.second);
            edges.push_back(e);
        }
        addEdgeLable(entry.first, entry.second, label);

        addEdge(entry.second, entry.first);
        addEdgeLable(entry.second, entry.first, label);
    }
    return edges;
}

std::unordered_set<std::string> Graph::getEdgeLabels(int v, int n)
{
    return edgeLables[v][n];
}

std::unordered_map<int, int> Graph::computingOrder()
{
    int maxDegree = 0;                  // 最大度数
    std::unordered_map<int, int> cores; // 存储每个节点的核心度
    std::vector<std::unordered_set<int>> orderedNodes(getNumberOfNodes());

    // 初始化orderedNodes
    for (int i = 0; i < getNumberOfNodes(); ++i)
    {
        orderedNodes[i] = std::unordered_set<int>();
    }

    // 将节点按度数分类
    std::unordered_map<int, int> degrees = getDegrees();
    for (const auto &entry : degrees)
    {
        orderedNodes[entry.second].insert(entry.first);
    }

    int node;             // 当前处理的节点
    int lowestDegree = 0; // 当前的最低度数
    int neighborDegree;   // 邻居的度数

    // 主循环
    while (lowestDegree < getNumberOfNodes())
    {
        if (orderedNodes[lowestDegree].empty())
        {
            ++lowestDegree; // 如果没有该度数的节点，增加度数
        }
        else
        {
            node = *orderedNodes[lowestDegree].begin();
            orderedNodes[lowestDegree].erase(node);
            cores[node] = lowestDegree; // 设置核心度
            degrees[node] = -1;         // 将节点度数设为-1，标记已处理

            // 更新所有邻居的度数
            for (int neighbor : getNeighbors(node))
            {
                neighborDegree = degrees[neighbor];
                if (neighborDegree > lowestDegree)
                {
                    orderedNodes[neighborDegree].erase(neighbor);
                    orderedNodes[neighborDegree - 1].insert(neighbor);
                    degrees[neighbor] = neighborDegree - 1;
                }
            }
        }
    }
    return cores;
}

int Graph::getNodeNum()
{
    return n;
}

std::unordered_set<std::string> Graph::checkLabelForm()
{
    std::unordered_set<std::string> labels;
    for(auto &entry : labelsForEdges){
        labels.insert(entry.first);
    }
    return labels;
}