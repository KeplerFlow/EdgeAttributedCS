#include "LCS.h"
LCS::LCS()
{
}

LCS::LCS(std::unordered_set<std::string> labels, Graph &graph)
{
    enumLabels(labels);
    // printSubLabels();
    build_in_mupt_index(graph);
    // opt(graph);
    //   printLcsIndex();
}

LCS::LCS(std::unordered_set<std::string> labels, Graph &graph, bool flag)
{
    enumLabels(labels);
    if (flag)
    {
        std::cout << "using the super optimization construction method." << std::endl;
        optWithSpend(graph);
    }
    else
    {
        std::cout << "using the normal construction method." << std::endl;
        build(graph);
    }
}

void LCS::incIndex(Graph &graph,int v,int u,std::string label){
    std::vector<std::unordered_set<std::string>> labels;
    for(std::unordered_set<std::string> set:subLabels){
        if(set.find(label)!=set.end()){
            labels.push_back(set);
        }
    }

    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForEdges = graph.getLabelsForEdges();

    //含特定标签的集合枚举
    for(std::unordered_set<std::string> set:labels){
        //仅用邻接表即可 邻接表根据标签重建;
        std::unordered_map<int, std::unordered_set<int>> adj;
        for(std::string label_temp:set){
            for(auto &entry:labelsForEdges[label_temp]){
                adj[entry.first].insert(entry.second);
                adj[entry.second].insert(entry.first);
            }
        }
        std::cout<<"增量维护索引,当前处理子标签集合: "<<std::endl;
        incIndex(adj,v,u,set,label);
    }
    std::cout<<"增量维护索引结束"<<std::endl;
}

void LCS::incIndex(std::unordered_map<int, std::unordered_set<int>> &adj,
                   int v, int u, std::unordered_set<std::string> labels, std::string label) {
    int root;
    if (vertexToK(u, labels) > vertexToK(v, labels)) {
        root = v;
    } else {
        root = u;
    }

    adj[u].insert(v);
    adj[v].insert(u);

    findSubCore(adj, root, labels);

    int k = vertexToK(root, labels);

    // 使用 std::set 来按 cd 从小到大排序
    std::set<std::pair<int, int>> sortedCd;

    // 初始化排序容器
    for (auto& p : cd) {
        sortedCd.insert({p.second, p.first}); // 插入 (cd值, 顶点) 对
    }

    bool flag = false;
    std::unordered_set<int> Vstar;
    
    while (!sortedCd.empty()) {
        auto pair = *sortedCd.begin(); // 获取最小的 cd 对
        sortedCd.erase(sortedCd.begin()); // 删除最小的 cd 对
        if(subcore.find(pair.first)==subcore.end()){
            continue;
        }
        if (flag) {
            Vstar.insert(pair.first);
            continue;
        }

        if (pair.first <= k) {
            for (int n : subcore[pair.second]) {
                if (cd[n] > pair.first) {
                    // 移除原来的（cd，顶点）对，并插入更新后的对
                    sortedCd.erase({cd[n], n});
                    cd[n]--;
                    sortedCd.insert({cd[n], n});
                }
            }
        } else {
            flag = true;
            Vstar.insert(pair.first);
        }
    }

    // 更新 LcsIndex 等部分
    for (int vertex : Vstar) {
        int core = vertexToK(vertex, labels);
        core++;
        for (int j = core; j > 0; j--)
        {
            if (LcsIndex[vertex][j].empty())
                {
                    LcsIndex[vertex][j].push_back(labels);
                    continue;
                }
            int i = 0;
            int flag = 2;
            while (i < LcsIndex[vertex][j].size())
            {
                std::unordered_set<std::string> labelss = LcsIndex[vertex][j][i];
                // 索引里的标签集是否为当前集合子集
                if (isSubsetOf(labels, labelss))
                {
                    flag = 0;
                    break;
                }
                // 当前集合是否为索引中标签集的子集
                if (isSubsetOf(labelss, labels))
                {
                    flag = 1;
                    LcsIndex[vertex][j].erase(LcsIndex[vertex][j].begin() + i);
                }
                i++;
            }
            if (flag != 0)
                LcsIndex[vertex][j].push_back(labels);
            else
                break;
        }
    }
}


void LCS::findSubCore(std::unordered_map<int, std::unordered_set<int>> &adj,int u,std::unordered_set<std::string> labels){
    //Init
    subcore.clear();
    std::queue<int> q;
    std::unordered_set<int> visited;
    int k = vertexToK(u,labels);

    q.push(u);
    visited.insert(u);

    //BFS
    while(!q.empty()){
        int v = q.front();
        q.pop();
        for(int n : adj[v]){
            int n_k = vertexToK(n,labels);
            if(n_k>=k){
                if(cd.find(v)==cd.end()){
                    cd[v]=0;
                }
                cd[v]++;
                if(n_k==k&&visited.find(n)==visited.end()){
                    q.push(n);
                    subcore[v].insert(n);
                    subcore[n].insert(v);
                    visited.insert(n);
                }
            }
        }
    }
}

int LCS::vertexToK(int u,std::unordered_set<std::string> labels){
    // Find the largest key
    int maxKey = std::max_element(LcsIndex[u].begin(), LcsIndex[u].end(), 
                                  [](const auto& lhs, const auto& rhs) {
                                      return lhs.first < rhs.first;
                                  })->first;

    // Check each key from maxKey down to 0
    for (int key = maxKey; key >= 0; --key) {
        auto it = LcsIndex[u].find(key);
        if (it != LcsIndex[u].end()) {
            // Check if any set in the vector is a subset of labels
            for (const auto& set : it->second) {
                if (isSubsetOf(labels, set)) {
                    return it->first;
                }
            }
        }
    }
    return -1;
}

void LCS::opt(Graph &graph)
{
    // 按照BFS顺序
    int size = 2;
    int place = 0;
    int right = 0, left = 0;
    int n = graph.getNodeNum() + 1;
    int edges = graph.getEdgeNum();
    int subnum = subLabels.size();
    // key is place

    std::unordered_map<int, std::shared_ptr<core::GLIST>> bfsCores(subnum);
    std::unordered_map<int, std::vector<std::vector<int>>> graphForCore(subnum);
    std::unordered_map<int, std::vector<int>> core(subnum);
    // key is size
    // std::unordered_map<int, std::unordered_set<std::unordered_set<std::string>>> labelSets;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForEdges = graph.getLabelsForEdges();
    std::cout << "begin to build the level; " << std::endl;
    int i = 0;
    for (auto &subset : subLabels)
    {
        if (subset.size() != 1)
            break;
        auto tempgraph = std::make_shared<Graph>(graph);
        tempgraph->filterEdgeLableAndRecover(subset);

        /**************/

        graphForCore[i].resize(n);
        core[i].resize(n);
        auto cm = std::make_shared<core::GLIST>(n);
        std::vector<std::pair<int, int>> edges = tempgraph->getEdges();

        for (int j = 0; j < edges.size(); ++j)
        {
            int v1 = edges.at(j).first;
            int v2 = edges.at(j).second;
            if (0 <= v1 && v1 < n && 0 <= v2 && v2 < n && v1 != v2)
            {
                graphForCore[i][v1].push_back(v2);
                graphForCore[i][v2].push_back(v1);
            }
        }
        cm->ComputeCore(graphForCore[i], true, core[i]);
        bfsCores[i] = cm;

        i++;
        right++;
    }

    for (left; left < right; left++)
    {
        std::vector<int> &core_left = core[left]; // 获取对应的vector引用，提高访问效率

        // 使用传统的for循环以获取索引和值
        for (size_t node = 0; node < core_left.size(); ++node)
        {
            int w = node; // 获取当前的值
            int k = core_left[node];

            // 检查从当前k值开始向下，是否存在包含当前labelSet的子集
            for (int i = k; i > 0; --i)
            {
                if (isSubSet(LcsIndex[w][i], subLabels[left]))
                {
                    break; // 一旦发现子集就停止检查
                }
                LcsIndex[w][i].push_back(subLabels[left]);
            }
        }
    }

    for (i; i < subLabels.size(); i++)
    {
        if (subLabels[i].size() == size)
        {
            // 找到上一层的子图
            int subIndex = 0;
            int min = INT_MAX;
            std::unordered_set<int> subIndexs;
            std::string targetLabel;
            // 遍历上层的图
            for (int j = 0; j < i; j++)
            {
                // 找到上一层符合要求的子图
                if (subLabels[j].size() == size - 1 && isSubsetOf(subLabels[i], subLabels[j]))
                {
                    // 找到多出的标签
                    for (auto &target : subLabels[i])
                    {
                        if (subLabels[j].find(target) == subLabels[j].end())
                        {
                            // 找到边最少的标签
                            if (min > labelsForEdges[target].size())
                            {
                                targetLabel = target;
                                subIndex = j;
                                min = labelsForEdges[target].size();
                            }
                            break;
                        }
                    }
                }
            }
            // 找到新加入的单个标签

            // 核维护得到新的图
            // std::vector<std::pair<int, int>> addedges = tempgraph->addEdgeSet(labelsForEdges[targetLabel], targetLabel);

            std::vector<std::vector<int>> tempEdges = graphForCore[subIndex];
            std::vector<std::pair<int, int>> edges;
            for (auto &entry : labelsForEdges[targetLabel])
            {
                int v1 = entry.first;
                int v2 = entry.second;
                if (entry.first != entry.second && find(tempEdges[entry.first].begin(), tempEdges[entry.first].end(), entry.second) == tempEdges[entry.first].end())
                {
                    if (v1 > v2)
                        std::swap(v1, v2);
                    const auto e = std::make_pair(v1, v2);
                    edges.push_back(e);
                    tempEdges[entry.first].push_back(entry.second);
                    tempEdges[entry.second].push_back(entry.first);
                }
            }

            // bfsGraphs[i] = tempgraph;

            auto cm = std::make_shared<core::GLIST>(bfsCores[subIndex]->exportGLISTData());

            graphForCore[i] = graphForCore[subIndex];
            core[i] = core[subIndex];

            // cm->ComputeCore(graphForCore[i], true, core[i]);

            for (int j = 0; j < edges.size(); j++)
            {
                cm->Insert(edges[j].first, edges[j].second, graphForCore[i], core[i]);
                // std::cout << "finish" << std::endl;
            }
            bfsCores[i] = cm;

            right++;
        }
        // 如果发现下一个标签集合隶属下一层BFS对象
        if (subLabels[i + 1].size() == size + 1)
        {
            for (left; left < right; left++)
            {
                std::vector<int> &core_left = core[left]; // 获取对应的vector引用，提高访问效率

                for (size_t node = 0; node < core_left.size(); ++node)
                {
                    int w = node; // 获取当前的值
                    int k = core_left[node];

                    // 检查从当前k值开始向下，是否存在包含当前labelSet的子集
                    for (int i = k; i > 0; --i)
                    {
                        if (isSubSet(LcsIndex[w][i], subLabels[left]))
                        {
                            break; // 一旦发现子集就停止检查
                        }
                        LcsIndex[w][i].push_back(subLabels[left]);
                    }
                }
            }

            std::cout << "next level : " << size << std::endl;
            int subIndex = 0;
            std::unordered_set<int> toErase;
            // 找到需要删除的标签集合索引
            for (subIndex = 0; subIndex < i; subIndex++)
            {
                if (subLabels[subIndex].size() == size - 1)
                    toErase.insert(subIndex);
            }
            // 删除相应的图信息
            for (int to : toErase)
            {
                // bfsGraphs.erase(to);
                //  bfsCoreness.erase(to);
                bfsCores.erase(to);
                core.erase(to);
                graphForCore.erase(to);
            }
            size++;
        }
    }
}

void LCS::opt_dfs(Graph &graph)
{
    // 一.按照DFS,进行初始化
    int size = 2;
    int place = 0;
    int right = 0, left = 0;
    int n = graph.getNodeNum() + 1;
    int edges = graph.getEdgeNum();
    int subnum = subLabels.size();
    // 为DFS选择上一层子图进行提前计算
    indexs_dfs.resize(subnum, 0);
    // 保存每个标签对应的边,用于核维护
    labelsForEdges_dfs = graph.getLabelsForEdges();
    std::unordered_set<std::string> labels = graph.checkLabelForm();

    // 预处理每个标签的边数
    std::unordered_map<std::string, int> edgeCount;
    for (auto &kv : labelsForEdges_dfs)
    {
        edgeCount[kv.first] = kv.second.size();
    }

    std::vector<std::pair<std::string, int>> sortedEdgeCount(edgeCount.begin(), edgeCount.end());

    // 将每个标签按照边数排序,使得核维护每次选择最少边的边标签加入
    std::sort(sortedEdgeCount.begin(), sortedEdgeCount.end(), [](const auto &a, const auto &b)
              {
                  return a.second < b.second; // 升序
              });

    // 计算第一层
    int temp_size = 2, temp_place = 0, temp_i = 0;
    for (auto &subset : subLabels)
    {
        if (subset.size() != 1)
        {
            break;
        }
        indexs_to_subset_dfs[temp_i]=subset;
        temp_i++;
    }
    // 从第二层开始,预计算从上一层哪个子标签进行核维护
    for (temp_i; temp_i < subLabels.size(); temp_i++)
    {
        if (subLabels[temp_i].size() == temp_size)
        {
            // 找到上一层的子图
            int subIndex = 0;
            int min = INT_MAX;
            std::unordered_set<int> subIndexs;
            std::string targetLabel;
            /*********** */
            for (const auto &kv : sortedEdgeCount)
            {
                if (subLabels[temp_i].find(kv.first) != subLabels[temp_i].end())
                {
                    std::unordered_set<std::string> sublas = subLabels[temp_i];
                    sublas.erase(kv.first);
                    for (int j = 0; j < temp_i; j++)
                    {
                        if (subLabels[j].size() == temp_size - 1 && isSubsetOf(sublas, subLabels[j]))
                        {
                            targetLabel = kv.first;
                            subIndex = j;
                            break;
                        }
                    }
                    break;
                }
            }
            /************* */
            indexs_dfs[subIndex]++;
            indexs_to_subset_dfs[temp_i]=subLabels[temp_i];
            treeIndexs_dfs[temp_i] = subIndex;
            treeLabelIndexs_dfs[temp_i] = targetLabel;
        }
        if (subLabels[temp_i + 1].size() == temp_size + 1)
        {
            temp_size++;
        }
    }

    // 先把第一层的算完,而后从第一层开始DFS
    std::cout << "begin to build the level in DFS; " << std::endl;
    int i = 0;
    for (auto &subset : subLabels)
    {
        if (subset.size() != 1)
            break;
        auto tempgraph = std::make_shared<Graph>(graph);
        tempgraph->filterEdgeLableAndRecover(subset);

        /**************/

        std::vector<std::vector<int>> graphForCore(n);
        std::vector<int> core(n);
        auto cm = std::make_shared<core::GLIST>(n);
        std::vector<std::pair<int, int>> edges = tempgraph->getEdges();

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

        // 将结果存储在LcsIndex中
        for (int node : core)
        {
            int k = core[node];

            for (int j = k; j > 0; j--)
            {
                if (LcsIndex[node][j].empty())
                {
                    LcsIndex[node][j].push_back(subset);
                    continue;
                }
                int i = 0;
                int flag = 2;
                while (i < LcsIndex[node][j].size())
                {
                    std::unordered_set<std::string> labelss = LcsIndex[node][j][i];
                    // 索引里的标签集是否为当前集合子集
                    if (isSubsetOf(subset, labelss))
                    {
                        flag = 0;
                        break;
                    }
                    // 当前集合是否为索引中标签集的子集
                    if (isSubsetOf(labelss, subset))
                    {
                        flag = 1;
                        LcsIndex[node][j].erase(LcsIndex[node][j].begin() + i);
                    }
                    i++;
                }
                if (flag != 0)
                    LcsIndex[node][j].push_back(subset);
                else
                    break;
            }
        }

        for (std::string label : labels)
        {
            if (subset.find(label) == subset.end())
            {
                if (indexs_dfs[i] != 0)
                {
                    subset.insert(label);
                    // DFS
                    subset.erase(label);
                }
            }
        }
        i++;
        right++;
    }
}

void LCS::dfs(Graph &graph, std::vector<std::vector<int>> graphForCore, std::vector<int> core, core::GLIST cm, std::unordered_set<std::string> &subset, int i)
{
    int subIndex = treeIndexs_dfs[i];
    std::string targetLabel = treeLabelIndexs_dfs[i];
    std::unordered_set<std::string> labels = graph.checkLabelForm();

    std::vector<std::vector<int>> tempEdges = graphForCore;
    std::vector<std::pair<int, int>> edges;
    for (auto &entry : labelsForEdges_dfs[targetLabel])
    {
        int v1 = entry.first;
        int v2 = entry.second;
        if (entry.first != entry.second && find(tempEdges[entry.first].begin(), tempEdges[entry.first].end(), entry.second) == tempEdges[entry.first].end())
        {
            if (v1 > v2)
                std::swap(v1, v2);
            const auto e = std::make_pair(v1, v2);
            edges.push_back(e);
            tempEdges[entry.first].push_back(entry.second);
            tempEdges[entry.second].push_back(entry.first);
        }
    }

    for (int j = 0; j < edges.size(); j++)
    {
        cm.Insert(edges[j].first, edges[j].second, graphForCore, core);
        // std::cout << "finish" << std::endl;
    }

    // 将结果存储在LcsIndex中
    for (int node : core)
    {
        int k = core[node];

        for (int j = k; j > 0; j--)
        {
            if (LcsIndex[node][j].empty())
            {
                LcsIndex[node][j].push_back(subset);
                continue;
            }
            int i = 0;
            int flag = 2;
            while (i < LcsIndex[node][j].size())
            {
                std::unordered_set<std::string> labelss = LcsIndex[node][j][i];
                // 索引里的标签集是否为当前集合子集
                if (isSubsetOf(subset, labelss))
                {
                    flag = 0;
                    break;
                }
                // 当前集合是否为索引中标签集的子集
                if (isSubsetOf(labelss, subset))
                {
                    flag = 1;
                    LcsIndex[node][j].erase(LcsIndex[node][j].begin() + i);
                }
                i++;
            }
            if (flag != 0)
                LcsIndex[node][j].push_back(subset);
            else
                break;
        }
    }

    for (std::string label : labels)
    {
        if (subset.find(label) == subset.end())
        {
            if (indexs_dfs[i] != 0)
            {
                subset.insert(label);
                // DFS
                subset.erase(label);
            }
        }
    }
}

void LCS::optWithSpend(Graph &graph)
{
    // 按照BFS顺序
    int size = 2;
    int place = 0;
    int right = 0, left = 0;
    int n = graph.getNodeNum() + 1;
    int edges = graph.getEdgeNum();
    int subnum = subLabels.size();
    // key is place

    std::unordered_map<int, std::shared_ptr<core::GLIST>> bfsCores(subnum);
    std::unordered_map<int, std::vector<std::vector<int>>> graphForCore(subnum);
    std::unordered_map<int, std::vector<int>> core(subnum);
    std::vector<int> indexs(subnum, 0);
    std::unordered_map<int, int> treeIndexs;
    std::unordered_map<int, std::string> treeLabelIndexs;
    // key is size
    // std::unordered_map<int, std::unordered_set<std::unordered_set<std::string>>> labelSets;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> labelsForEdges = graph.getLabelsForEdges();

    // 预处理每个标签的边数
    std::unordered_map<std::string, int> edgeCount;
    for (auto &kv : labelsForEdges)
    {
        edgeCount[kv.first] = kv.second.size();
    }

    // 将哈希表转换为可排序的vector
    std::vector<std::pair<std::string, int>> sortedEdgeCount(edgeCount.begin(), edgeCount.end());

    // 按值排序，使用自定义的比较函数
    std::sort(sortedEdgeCount.begin(), sortedEdgeCount.end(), [](const auto &a, const auto &b)
              {
                  return a.second < b.second; // 升序
              });

    // 记录要存储的上一层信息
    int temp_size = 2, temp_place = 0, temp_i = 0;
    for (auto &subset : subLabels)
    {
        if (subset.size() != 1)
        {
            break;
        }
        temp_i++;
    }
    for (temp_i; temp_i < subLabels.size(); temp_i++)
    {
        if (subLabels[temp_i].size() == temp_size)
        {
            // 找到上一层的子图
            int subIndex = 0;
            int min = INT_MAX;
            std::unordered_set<int> subIndexs;
            std::string targetLabel;
            /*********** */
            for (const auto &kv : sortedEdgeCount)
            {
                if (subLabels[temp_i].find(kv.first) != subLabels[temp_i].end())
                {
                    std::unordered_set<std::string> sublas = subLabels[temp_i];
                    sublas.erase(kv.first);
                    for (int j = 0; j < temp_i; j++)
                    {
                        if (subLabels[j].size() == temp_size - 1 && isSubsetOf(sublas, subLabels[j]))
                        {
                            targetLabel = kv.first;
                            subIndex = j;
                            break;
                        }
                    }
                    break;
                }
            }
            /************* */
            indexs[subIndex]++;
            treeIndexs[temp_i] = subIndex;
            treeLabelIndexs[temp_i] = targetLabel;
        }
        if (subLabels[temp_i + 1].size() == temp_size + 1)
        {
            temp_size++;
        }
    }

    std::cout << "begin to build the level; " << std::endl;
    int i = 0;
    for (auto &subset : subLabels)
    {
        if (subset.size() != 1)
            break;
        auto tempgraph = std::make_shared<Graph>(graph);
        tempgraph->filterEdgeLableAndRecover(subset);

        /**************/

        graphForCore[i].resize(n);
        core[i].resize(n);
        auto cm = std::make_shared<core::GLIST>(n);
        std::vector<std::pair<int, int>> edges = tempgraph->getEdges();

        for (int j = 0; j < edges.size(); ++j)
        {
            int v1 = edges.at(j).first;
            int v2 = edges.at(j).second;
            if (0 <= v1 && v1 < n && 0 <= v2 && v2 < n && v1 != v2)
            {
                graphForCore[i][v1].push_back(v2);
                graphForCore[i][v2].push_back(v1);
            }
        }
        cm->ComputeCore(graphForCore[i], true, core[i]);
        if (indexs[i] != 0)
            bfsCores[i] = cm;

        i++;
        right++;
    }

    for (left; left < right; left++)
    {
        std::vector<int> &core_left = core[left]; // 获取对应的vector引用，提高访问效率

        for (size_t node = 0; node < core_left.size(); ++node)
        {
            int w = node; // 获取当前的值
            int k = core_left[node];

            // 检查从当前k值开始向下，是否存在包含当前labelSet的子集
            for (int i = k; i > 0; --i)
            {
                if (isSubSet(LcsIndex[w][i], subLabels[left]))
                {
                    break; // 一旦发现子集就停止检查
                }
                LcsIndex[w][i].push_back(subLabels[left]);
            }
        }
    }

    for (i; i < subLabels.size(); i++)
    {
        if (subLabels[i].size() == size)
        {
            // 找到上一层的子图
            int subIndex = treeIndexs[i];
            std::string targetLabel = treeLabelIndexs[i];

            std::vector<std::vector<int>> tempEdges = graphForCore[subIndex];
            std::vector<std::pair<int, int>> edges;
            for (auto &entry : labelsForEdges[targetLabel])
            {
                int v1 = entry.first;
                int v2 = entry.second;
                if (entry.first != entry.second && find(tempEdges[entry.first].begin(), tempEdges[entry.first].end(), entry.second) == tempEdges[entry.first].end())
                {
                    if (v1 > v2)
                        std::swap(v1, v2);
                    const auto e = std::make_pair(v1, v2);
                    edges.push_back(e);
                    tempEdges[entry.first].push_back(entry.second);
                    tempEdges[entry.second].push_back(entry.first);
                }
            }

            auto cm = std::make_shared<core::GLIST>(bfsCores[subIndex]->exportGLISTData());
            graphForCore[i] = graphForCore[subIndex];
            core[i] = core[subIndex];

            for (int j = 0; j < edges.size(); j++)
            {
                cm->Insert(edges[j].first, edges[j].second, graphForCore[i], core[i]);
                // std::cout << "finish" << std::endl;
            }

            if (indexs[i] != 0)
            {
                bfsCores[i] = cm;
            }

            right++;
        }
        // 如果发现下一个标签集合隶属下一层BFS对象
        if (subLabels[i + 1].size() == size + 1)
        {
            for (left; left < right; left++)
            {
                std::vector<int> &core_left = core[left]; // 获取对应的vector引用，提高访问效率

                for (size_t node = 0; node < core_left.size(); ++node)
                {
                    int w = node; // 获取当前的值
                    int k = core_left[node];

                    // 检查从当前k值开始向下，是否存在包含当前labelSet的子集
                    for (int i = k; i > 0; --i)
                    {
                        if (isSubSet(LcsIndex[w][i], subLabels[left]))
                        {
                            break; // 一旦发现子集就停止检查
                        }
                        LcsIndex[w][i].push_back(subLabels[left]);
                    }
                }
            }

            std::cout << "next level : " << size << std::endl;
            int subIndex = 0;
            std::unordered_set<int> toErase;
            // 找到需要删除的标签集合索引
            for (subIndex = 0; subIndex < i; subIndex++)
            {
                if (subLabels[subIndex].size() == size - 1)
                    toErase.insert(subIndex);
            }
            // 删除相应的图信息
            for (int to : toErase)
            {
                // bfsGraphs.erase(to);
                //  bfsCoreness.erase(to);
                bfsCores.erase(to);
                core.erase(to);
                graphForCore.erase(to);
            }
            size++;
        }
    }
}

void LCS::build(Graph &graph)
{
    for (auto &labelSet : normalSubLabels)
    {
        // 存储被删除的边
        Graph tempgraph(graph);
        // tempgraph.printGraphInfo(tempgraph);

        tempgraph.filterEdgeLableAndRecover(labelSet);
        tempgraph.computeDegrees();
        tempgraph.computeMinimumDegree();

        std::unordered_map<int, int> coreness = coreDecomposition(labelSet, tempgraph);

        // 将结果存储在LcsIndex中
        for (const auto &pair : coreness)
        {
            int node = pair.first;
            int k = pair.second;

            for (int j = k; j > 0; j--)
            {
                if (LcsIndex[node][j].empty())
                {
                    LcsIndex[node][j].push_back(labelSet);
                    continue;
                }
                int i = 0;
                int flag = 2;
                while (i < LcsIndex[node][j].size())
                {
                    std::unordered_set<std::string> labels = LcsIndex[node][j][i];
                    // 索引里的标签集是否为当前集合子集
                    if (isSubsetOf(labelSet, labels))
                    {
                        flag = 0;
                        break;
                    }
                    // 当前集合是否为索引中标签集的子集
                    if (isSubsetOf(labels, labelSet))
                    {
                        flag = 1;
                        LcsIndex[node][j].erase(LcsIndex[node][j].begin() + i);
                    }
                    i++;
                }
                if (flag != 0)
                    LcsIndex[node][j].push_back(labelSet);
                else
                    break;
            }
        }
    }
}

void LCS::build_in_mupt_index(Graph &graph)
{
    std::chrono::milliseconds totalDuration(0);

    for (int i=0;i<normalSubLabels .size();i++)
    {
        auto lcsStart = std::chrono::high_resolution_clock::now();
        
        std::unordered_set<std::string> labelSet = normalSubLabels [i];
        // 存储被删除的边
        
        Graph temp(graph);
        temp.filterEdgeLableAndRecover(labelSet);
        temp.computeDegrees();
        temp.computeMinimumDegree();
        std::unordered_map<int, std::unordered_set<int>> cores;
        std::unordered_map<int, int> core = CoreGroup::coreGroupsAlgorithmInSet(temp,cores);
            
        auto mid = std::chrono::high_resolution_clock::now();
        mupt_index[i] = cores;
        std::chrono::milliseconds loopDuration = std::chrono::duration_cast<std::chrono::milliseconds>(mid - lcsStart);
        
        // 累加时间
        totalDuration += loopDuration;
    }

    std::cout << "the basic LCS Index built in " << totalDuration.count() << " ms" << std::endl;

}

void LCS::transIndex()
{
    for (int i=0;i<normalSubLabels.size();i++){
        for (auto &pair: mupt_index[i]){
            for (auto &v:pair.second){
                mupt_index_map[i].emplace(pair.first, v);
            }
        }
    }
}

std::unordered_map<int,std::unordered_set<int>> LCS::mupt_index_search(std::unordered_set<std::string> labels){
    for (int i=0;i<normalSubLabels.size();i++){
        if(normalSubLabels[i].size()==labels.size()&&isSubsetOf(labels,normalSubLabels[i])){
            return mupt_index[i];
        }
    }
    return {};
}

std::unordered_set<int> LCS::mupt_index_search_map(std::unordered_set<std::string> labels,int k){
    int index=-1;
    std::unordered_set<int> solution;
    for (int i=0;i<normalSubLabels.size();i++){
        if(normalSubLabels[i].size()==labels.size()&&isSubsetOf(labels,normalSubLabels[i])){
            index = i;
            break;
        }
    }
    for (auto it = mupt_index_map[index].begin(); it != mupt_index_map[index].end(); ++it) {
        if (it->first < k) {
            break; // 当 key 小于 k 时停止遍历
        }
        solution.insert(it->second); // 插入当前 vertex ID 到 solution
    }
    return solution;
}

// 根据边标签集截断,得到每个顶点对应的coreness
std::unordered_map<int, int> LCS::coreDecomposition(std::unordered_set<std::string> labels, Graph &graph)
{
    return CoreGroup::coreGroupsAlgorithm(graph);
}

bool LCS::isExist(int vertex, std::unordered_set<std::string> &labels, int k)
{
    return isSubSet(LcsIndex[vertex][k], labels);
}

// 前者的每个集合中是否存在后者的子集
bool LCS::isSubSet(std::vector<std::unordered_set<std::string>> &lcs, std::unordered_set<std::string> &labelSet)
{
    // 遍历 lcs 中的每个集合
    for (const auto &set : lcs)
    {
        // 即set是否是labelSet的子集?
        if (isSubsetOf(labelSet, set))
        {
            return true; // 找到一个子集
        }
    }
    return false; // 没有找到任何子集
}

// 后者是否为前者子集
bool LCS::isSubsetOf(const std::unordered_set<std::string> &subset, const std::unordered_set<std::string> &set)
{
    return std::all_of(set.begin(), set.end(), [&](const std::string &item)
                       { return subset.find(item) != subset.end(); });
}

// 枚举边标签集得到子集
void LCS::enumLabels(std::unordered_set<std::string> labels)
{
    std::vector<std::string> labelList(labels.begin(), labels.end());
    int n = labelList.size();
    std::vector<std::unordered_set<std::string>> tempSubLabels;

    // 从 1 到 2^n - 1 遍历所有可能的子集
    for (int i = 1; i < (1 << n); ++i)
    {
        std::unordered_set<std::string> subset;
        for (int j = 0; j < n; ++j)
        {
            // 检查第 j 位是否被设置
            if (i & (1 << j))
            {
                subset.insert(labelList[j]);
            }
        }
        if (!subset.empty())
        {
            tempSubLabels.push_back(subset);
            bfsLabels[subset.size()].push_back(subset);
        }
    }
    normalSubLabels = tempSubLabels;

    // 根据集合大小进行排序
    std::sort(tempSubLabels.begin(), tempSubLabels.end(), [](const std::unordered_set<std::string> &a, const std::unordered_set<std::string> &b)
              { return a.size() < b.size(); });

    // 存储排序后的结果
    subLabels = std::move(tempSubLabels);
}

void LCS::printSubLabels()
{
    for (auto &subset : subLabels)
    {
        for (auto &elem : subset)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void LCS::printLcsIndex()
{
    for (const auto &vertex : LcsIndex)
    {
        std::cout << "Vertex ID: " << vertex.first << std::endl;
        for (const auto &k : vertex.second)
        {
            std::cout << "  K value: " << k.first << std::endl;
            int setIndex = 0;
            for (const auto &edgeSet : k.second)
            {
                std::cout << "    Set " << ++setIndex << ": ";
                for (const std::string &label : edgeSet)
                {
                    std::cout << label << " ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << "-----------------------------------------------" << std::endl;
    }
}
// 计算 LcsIndex 占用的内存
size_t LCS::calculateMemoryUsage()
{
    size_t total_memory = 0;
    double k_sum = 0;
    double sum_of_set = 0;
    // 遍历外部的 unordered_map
    for (const auto &outer_pair : LcsIndex)
    {
        //total_memory += sizeof(outer_pair.first);  // 顶点id大小
        //total_memory += sizeof(outer_pair.second); // 内部unordered_map大小
        k_sum+=outer_pair.second.size();
        // 遍历内部的 unordered_map
        for (const auto &inner_pair : outer_pair.second)
        {
            //total_memory += sizeof(inner_pair.first);  // k值大小
            //total_memory += sizeof(inner_pair.second); // vector大小
            sum_of_set+=inner_pair.second.size();
            // 遍历 vector
            for (const auto &label_set : inner_pair.second)
            {
                //total_memory += sizeof(label_set); // unordered_set大小

                // 遍历 unordered_set
                for (const auto &label : label_set)
                {
                    total_memory += sizeof(char);    // std::string 对象本身的大小
                    //total_memory += label.capacity(); // 字符串内容占用的大小
                }
            }
        }
    }
    double l = sum_of_set/(k_sum);
    std::cout<<"每个单元格的平均集合数量为 "<<l<<std::endl;
    return total_memory;
}
void LCS::serialize(std::ofstream &ofs) const
{
    // 序列化 LcsIndex
    size_t outer_map_size = LcsIndex.size();
    ofs.write(reinterpret_cast<const char *>(&outer_map_size), sizeof(outer_map_size));
    for (const auto &outer_pair : LcsIndex)
    {
        ofs.write(reinterpret_cast<const char *>(&outer_pair.first), sizeof(outer_pair.first)); // 顶点id
        size_t inner_map_size = outer_pair.second.size();
        ofs.write(reinterpret_cast<const char *>(&inner_map_size), sizeof(inner_map_size));
        for (const auto &inner_pair : outer_pair.second)
        {
            ofs.write(reinterpret_cast<const char *>(&inner_pair.first), sizeof(inner_pair.first)); // k值
            size_t vector_size = inner_pair.second.size();
            ofs.write(reinterpret_cast<const char *>(&vector_size), sizeof(vector_size));
            for (const auto &label_set : inner_pair.second)
            {
                size_t set_size = label_set.size();
                ofs.write(reinterpret_cast<const char *>(&set_size), sizeof(set_size));
                for (const auto &label : label_set)
                {
                    size_t label_size = label.size();
                    ofs.write(reinterpret_cast<const char *>(&label_size), sizeof(label_size));
                    ofs.write(label.c_str(), label_size); // 写入标签字符串
                }
            }
        }
    }

   
}

void LCS::serialize_mupt_index(std::ofstream &ofs) const {
    // 序列化 mupt_index
    size_t outer_map_size = mupt_index.size();
    ofs.write(reinterpret_cast<const char *>(&outer_map_size), sizeof(outer_map_size));

    for (const auto &outer_pair : mupt_index) {
        ofs.write(reinterpret_cast<const char *>(&outer_pair.first), sizeof(outer_pair.first)); // 写入外层 key (int)

        size_t inner_map_size = outer_pair.second.size();
        ofs.write(reinterpret_cast<const char *>(&inner_map_size), sizeof(inner_map_size)); // 写入内层 map 大小

        for (const auto &inner_pair : outer_pair.second) {
            ofs.write(reinterpret_cast<const char *>(&inner_pair.first), sizeof(inner_pair.first)); // 写入内层 key (int)

            size_t set_size = inner_pair.second.size();
            ofs.write(reinterpret_cast<const char *>(&set_size), sizeof(set_size)); // 写入 set 大小

            for (const auto &value : inner_pair.second) {
                ofs.write(reinterpret_cast<const char *>(&value), sizeof(value)); // 写入 set 中的值 (int)
            }
        }
    }
}

void LCS::deserialize_mupt_index(std::ifstream &ifs) {
    // 反序列化 mupt_index
    size_t outer_map_size;
    ifs.read(reinterpret_cast<char *>(&outer_map_size), sizeof(outer_map_size));

    for (size_t i = 0; i < outer_map_size; ++i) {
        int outer_key;
        ifs.read(reinterpret_cast<char *>(&outer_key), sizeof(outer_key)); // 读取外层 key

        size_t inner_map_size;
        ifs.read(reinterpret_cast<char *>(&inner_map_size), sizeof(inner_map_size)); // 读取内层 map 大小

        std::unordered_map<int, std::unordered_set<int>> inner_map;
        for (size_t j = 0; j < inner_map_size; ++j) {
            int inner_key;
            ifs.read(reinterpret_cast<char *>(&inner_key), sizeof(inner_key)); // 读取内层 key

            size_t set_size;
            ifs.read(reinterpret_cast<char *>(&set_size), sizeof(set_size)); // 读取 set 大小

            std::unordered_set<int> value_set;
            for (size_t k = 0; k < set_size; ++k) {
                int value;
                ifs.read(reinterpret_cast<char *>(&value), sizeof(value)); // 读取 set 中的值
                value_set.insert(value);
            }
            inner_map[inner_key] = value_set;
        }
        mupt_index[outer_key] = inner_map;
    }
}

void LCS::deserialize(std::ifstream &ifs)
{
    // 反序列化 LcsIndex
    size_t outer_map_size;
    ifs.read(reinterpret_cast<char *>(&outer_map_size), sizeof(outer_map_size));
    for (size_t i = 0; i < outer_map_size; ++i)
    {
        int vertex_id;
        ifs.read(reinterpret_cast<char *>(&vertex_id), sizeof(vertex_id));

        size_t inner_map_size;
        ifs.read(reinterpret_cast<char *>(&inner_map_size), sizeof(inner_map_size));

        std::unordered_map<int, std::vector<std::unordered_set<std::string>>> inner_map;
        for (size_t j = 0; j < inner_map_size; ++j)
        {
            int k_value;
            ifs.read(reinterpret_cast<char *>(&k_value), sizeof(k_value));

            size_t vector_size;
            ifs.read(reinterpret_cast<char *>(&vector_size), sizeof(vector_size));

            std::vector<std::unordered_set<std::string>> label_sets;
            for (size_t k = 0; k < vector_size; ++k)
            {
                size_t set_size;
                ifs.read(reinterpret_cast<char *>(&set_size), sizeof(set_size));

                std::unordered_set<std::string> label_set;
                for (size_t l = 0; l < set_size; ++l)
                {
                    size_t label_size;
                    ifs.read(reinterpret_cast<char *>(&label_size), sizeof(label_size));

                    std::string label(label_size, '\0');
                    ifs.read(&label[0], label_size);
                    label_set.insert(label);
                }
                label_sets.push_back(label_set);
            }
            inner_map[k_value] = label_sets;
        }
        LcsIndex[vertex_id] = inner_map;
    }
}
