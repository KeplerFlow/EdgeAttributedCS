#include "hierarchy.h"

Hierarchy::Hierarchy(Graph &graph)
{
	nodeSum = 0;
	kmax = 0;
	// core decompositon
	vertexSum = coreDecomposition(graph);
	build(graph);
	for (auto &pair : coreIndex)
	{
		// std::cout << "the vertex is " << pair.first << " the coreness is " << pair.second << std::endl;
		coreness[pair.second].insert(pair.first);
		if (pair.second > kmax)
		{
			kmax = pair.second;
		}
	}
	// printData();
	//  std::unordered_set<int> checksum;
	//  for (auto &pair : nodes)
	//  {
	//  	checksum.insert(pair.second.begin(), pair.second.end());
	//  }
	//  std::unordered_set<int> toFind;
	//  std::unordered_set<int> graphNodes = graph.getNodes();
	//  for (const auto &node : graphNodes)
	//  {
	//  	if (checksum.find(node) == checksum.end())
	//  	{
	//  		toFind.insert(node);
	//  		std::cout << "the node " << node << " is not in the hierarchy and the coreness is " << coreIndex[node] << std::endl;
	//  	}
	//  }
	//  std::cout << "the checksum is " << checksum.size() << std::endl;
	//  std::cout << "the node size is " << graph.getNodes().size() << std::endl;
	//  std::cout << "the to find size is " << toFind.size() << std::endl;
	//  print vertexToNode
	//  for (auto &pair : vertexToNode)
	//  {
	//  	std::cout << "the vertex is " << pair.first << " the node is " << pair.second << std::endl;
	//  }
}

int Hierarchy::coreDecomposition(Graph &graph)
{
	coreIndex = CoreGroup::coreGroupsAlgorithm(graph);
	for (auto &pair : coreIndex)
	{
		// std::cout << "the vertex is " << pair.first << " the coreness is " << pair.second << std::endl;
		coreness[pair.second].insert(pair.first);
		if (pair.second > kmax)
		{
			kmax = pair.second;
		}
	}
	std::cout << "the kmax is " << kmax << std::endl;
	return coreIndex.size();
}

void Hierarchy::build(Graph &graph)
{
	// Initialize an empty root
	nodes[nodeSum] = {}; // Let's assume root is at index 0
	nodeSum++;

	std::unordered_set<int> toErase;
	for (int vertex : coreness[1])
	{
		toErase.insert(vertex);
		vertexToNode[vertex] = nodeSum;
	}
	for (int vertex : toErase)
	{
		coreness[1].erase(vertex);
	}

	nodes[nodeSum] = toErase;
	parents[nodeSum] = nodeSum - 1;
	children[nodeSum - 1].insert(nodeSum);
	// std::unordered_set<std::string> labels;
	// for(int v:nodes[nodeSum]){
	// 	for(int neighbor : graph.getNeighbors(v)){
	// 		if(nodes[nodeSum-1].find(neighbor)!=nodes[nodeSum-1].end()){
	// 			for(std::string label:graph.getEdgeLabels(v,neighbor))
	// 			{
	// 				labels.insert(label);
	// 			}
	// 		}
	// 	}
	// }
	// edgeLabel[nodeSum][nodeSum-1]=labels;
	// edgeLabel[nodeSum-1][nodeSum]=labels;

	nodeSum++;
	std::unordered_map<int, bool> visited;
	int minCoreness = 2;
	// 从2-core开始才会分裂
	DFS(graph, coreness, minCoreness, 1);

	// 继而建立节点间标签集合信息
	// for (int i = 1; i < nodeSum; i++)
	// {
	// 	std::unordered_set<std::string> labels;
	// 	int par = parents[i];
	// 	for (int v : nodes[i])
	// 	{
	// 		for (int neighbor : graph.getNeighbors(v))
	// 		{
	// 			if (nodes[par].find(neighbor) != nodes[par].end())
	// 			{
	// 				for (std::string label : graph.getEdgeLabels(v, neighbor))
	// 				{
	// 					labels.insert(label);
	// 				}
	// 			}
	// 		}
	// 	}
	// 	edgeLabel[i][par] = labels;
	// 	edgeLabel[par][i] = labels;
	// 	std::cout << "the size of this " << labels.size() << std::endl;
	// }

	for (int i = 1; i < nodeSum; i++)
	{
		for (int child : children[i])
		{
			std::unordered_set<std::string> labels;
			std::unordered_set<int> ver;

			std::queue<int> q;
			q.push(child);
			while (!q.empty())
			{
				int temp = q.front();
				for (int v : nodes[temp])
				{
					ver.insert(v);
				}
				q.pop();
				for (int ch : children[temp])
				{
					q.push(ch);
				}
			}

			for (int v : ver)
			{
				for (int neighbor : graph.getNeighbors(v))
				{
					if (nodes[child].find(neighbor) != nodes[child].end())
					{
						for (std::string label : graph.getEdgeLabels(v, neighbor))
						{
							labels.insert(label);
						}
					}
				}
			}
			edgeLabel[i][child] = labels;
			edgeLabel[child][i] = labels;
		}
	}
}

// 根据传入的高core顶点,选择
void Hierarchy::DFS(Graph &graph, std::unordered_map<int, std::unordered_set<int>> &vertexs, int minCoreness, int nodeId)
{
	assert(vertexs.size() > 0 && "vertexs map should not be empty");
	assert(vertexSum >= 0 && "vertexSum should be greater than 0");

	std::unordered_map<int, std::unordered_set<int>> localVertexs = vertexs; // 使用本地副本
	while (!localVertexs.empty() && vertexSum >= 0)
	{
		while (localVertexs[minCoreness].empty() && minCoreness < kmax)
		{
			minCoreness++;
		}
		if (localVertexs[minCoreness].empty())
			break;

		int vertex = *localVertexs[minCoreness].begin();
		if (coreness[minCoreness].find(vertex) == coreness[minCoreness].end())
		{
			localVertexs[minCoreness].erase(vertex);
			if (localVertexs[minCoreness].empty())
				localVertexs.erase(minCoreness);
			continue;
		}
		int curNodeId = 0;
		std::unordered_map<int, std::unordered_set<int>> newVertexs = BFS(graph, vertex, localVertexs, minCoreness, nodeId, curNodeId);

		if (!newVertexs.empty())
			// std::cout << "new vertexs is not empty,Keeping DFS" << std::endl;
			if (!newVertexs.empty())
				DFS(graph, newVertexs, minCoreness + 1, curNodeId);
	}
}

// nodeid为该节点的父节点id
// curNodeId为该节点的id
std::unordered_map<int, std::unordered_set<int>> Hierarchy::BFS(Graph &graph, int startVertex, std::unordered_map<int, std::unordered_set<int>> &vertexs, int minCoreness, int nodeId, int &curNodeId)
{

	// std::cout << "start vertex is " << startVertex << std::endl;
	std::unordered_map<int, std::unordered_set<int>> component;
	std::unordered_set<int> toErase;

	std::unordered_set<int> visited;
	std::queue<int> queue;
	queue.push(startVertex);
	visited.insert(startVertex);
	toErase.insert(startVertex);
	vertexToNode[startVertex] = nodeSum;
	while (!queue.empty())
	{
		int current = queue.front();
		queue.pop();

		for (int neighbor : graph.getNeighbors(current))
		{
			if (coreIndex[neighbor] >= minCoreness)
			{
				if (vertexs[coreIndex[neighbor]].find(neighbor) != vertexs[coreIndex[neighbor]].end())
				{
					if (coreIndex[neighbor] == minCoreness && visited.find(neighbor) == visited.end())
					{
						queue.push(neighbor);
						visited.insert(neighbor);

						vertexToNode[neighbor] = nodeSum;
						toErase.insert(neighbor);
					}
					else if (coreIndex[neighbor] > minCoreness && visited.find(neighbor) == visited.end())
					{
						queue.push(neighbor);
						visited.insert(neighbor);

						component[coreIndex[neighbor]].insert(neighbor);
					}
				}
			}
		}
	}
	for (int item : toErase)
	{
		// std::cout << "at the " << minCoreness << "-core , erase vertex " << item << std::endl;

		coreness[coreIndex[item]].erase(item);

		if (coreness[coreIndex[item]].empty())
		{
			coreness.erase(coreIndex[item]);
			// std::cout << "erase the " << coreIndex[item] << "-core" << std::endl;
		}
		vertexSum--;
	}
	// std::cout << "left vertex sum is " << vertexSum << std::endl;
	//  Node node(nodeId, minCoreness, res);
	nodes[nodeSum] = toErase;
	parents[nodeSum] = nodeId;
	nodeCoreness[nodeSum] = minCoreness;
	// std::cout << "the current node is " << nodeSum << " and the parent node is " << nodeId << std::endl;
	children[nodeId].insert(nodeSum);
	// std::cout << "the current node is " << nodeId << " and the children node is " << nodeSum << std::endl;
	curNodeId = nodeSum;

	// std::unordered_set<std::string> labels;
	// for(int v:nodes[nodeSum]){
	// 	for(int neighbor : graph.getNeighbors(v)){
	// 		if(nodes[nodeId].find(neighbor)!=nodes[nodeId].end()){
	// 			for(std::string label:graph.getEdgeLabels(v,neighbor))
	// 			{
	// 				labels.insert(label);
	// 			}
	// 		}
	// 	}
	// }
	// edgeLabel[nodeSum][nodeId]=labels;
	// edgeLabel[nodeId][nodeSum]=labels;

	nodeSum++;
	return component;
}
std::unordered_set<int> Hierarchy::query(int queryNodes, int k)
{
	std::unordered_set<int> result;
	if (coreIndex[queryNodes] < k)
	{
		// std::cout << "the query node is not in the k-core" << std::endl;
		return {};
	}
	int targetNodeId = vertexToNode[queryNodes];
	// std::cout << "the query node is " << queryNodes << " the targetNodeId is " << targetNodeId << std::endl;
	// std::cout << "the target node's coreness is " << nodeCoreness[targetNodeId] << std::endl;
	while (nodeCoreness[targetNodeId] >= k)
	{
		// 如果父节点的核心度小于k，那么就是当前节点id就是子树的根
		if (nodeCoreness[parents[targetNodeId]] < k)
		{
			// std::cout << "break,the target node is " << targetNodeId << std::endl;
			break;
		} // 如果父节点的核心度等于k，那么就是父节点的id就是子树的根
		else if (nodeCoreness[parents[targetNodeId]] = k)
		{
			// std::cout << "the target node is " << targetNodeId << std::endl;
			targetNodeId = parents[targetNodeId];
		} // 如果父节点的核心度大于k，那么先到父节点
		else
		{
			// std::cout << "the target node is " << targetNodeId << std::endl;
			targetNodeId = parents[parents[targetNodeId]];
		}
	}
	// 从targetNodeId开始BFS遍历子树
	std::queue<int> queue;
	queue.push(targetNodeId);
	while (!queue.empty())
	{
		int current = queue.front();
		for (int vertex : nodes[current])
		{
			result.insert(vertex);
		}
		// std::cout << "the current node is " << current << std::endl;
		queue.pop();
		for (int child : children[current])
		{
			queue.push(child);
		}
	}
	return result;
}

std::unordered_set<int> Hierarchy::bl(int querynode, int k, Graph &tempgraph, std::unordered_set<std::string> &labels)
{
	std::unordered_set<int> result;
	if (coreIndex[querynode] < k)
	{
		// std::cout << "the query node is not in the k-core" << std::endl;
		return {};
	}
	int targetNodeId = vertexToNode[querynode];

	while (nodeCoreness[targetNodeId] >= k)
	{
		int parentNodeId = parents[targetNodeId];
		if (nodeCoreness[parentNodeId] <= k)
		{
			if (nodeCoreness[parentNodeId] < k)
			{
				break;
			}
			if (isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
			{
				targetNodeId = parentNodeId;
			}
			else
			{
				break;
			}
			// targetNodeId = parentNodeId;
		}
		else
		{
			if (isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
			{
				targetNodeId = parentNodeId;
			}
			else
			{
				break;
			}
		}
	}
	// 从targetNodeId开始BFS遍历子树
	std::queue<int> queue;
	queue.push(targetNodeId);
	while (!queue.empty())
	{
		int current = queue.front();
		for (int vertex : nodes[current])
		{
			result.insert(vertex);
		}
		// std::cout << "the current node is " << current << std::endl;
		queue.pop();
		for (int child : children[current])
		{
			queue.push(child);
		}
	}
	std::cout << "the CL-tree size is " << result.size() << std::endl;

	// 得到度数-节点集合 节点-度数集合 节点集合
	std::unordered_map<int, std::unordered_set<int>> largestComponent;
	int i = 0;
	for (std::string label : labels)
	{
		std::unordered_set<std::string> s;
		s.insert(label);
		Graph temp(tempgraph);
		temp.filterEdgeLableAndRecover(s);
		std::vector<std::unordered_set<int>> degreeVec = temp.getOrderedNodes();
		std::unordered_map<int, int> nodeDegrees = temp.getDegrees();
		std::unordered_set<int> community = result;
		// 删除节点直到最小度数达到或超过k
		// std::cout << "度数重建完毕 " << label << std::endl;
		int flag = false;

		int minDegree = 0;

		int nodeToDelete = -1;

		// 找到第一个待删除的节点
		for (int j = 0; j < k; j++)
		{
			if (!degreeVec[j].empty())
			{
				flag = true;
				minDegree = j;
				nodeToDelete = *degreeVec[j].begin();
				break;
			}
		}
		// std::cout << "开启顶点选择完毕 " << label << std::endl;
		while (flag)
		{
			// std::cout << community.size() << std::endl;
			// std::cout<<"the node to delete is "<<nodeToDelete<<std::endl;
			community.erase(nodeToDelete);
			degreeVec[nodeDegrees[nodeToDelete]].erase(nodeToDelete);
			for (int neighbor : temp.getNeighbors(nodeToDelete))
			{
				if (community.find(neighbor) != community.end())
				{ // 仅更新仍在社区中的节点
					degreeVec[nodeDegrees[neighbor]].erase(neighbor);
					if (nodeDegrees[neighbor] > 0)
					{
						nodeDegrees[neighbor]--;
						degreeVec[nodeDegrees[neighbor]].insert(neighbor);
					}
				}
			}
			flag = false;
			// 找到第一个待删除的节点
			for (int j = 0; j < k; j++)
			{
				if (!degreeVec[j].empty())
				{
					flag = true;
					minDegree = j;
					nodeToDelete = *degreeVec[j].begin();
					break;
				}
			}
		} // 从查询节点开始检索最大连通分量
		largestComponent[i] = community;

		i++;
	}

	// std::cout << "开始合并" << std::endl;

	std::unordered_set<int> intersection;

	for (const auto &elem : largestComponent[0])
	{
		if (largestComponent[1].find(elem) != largestComponent[1].end())
		{
			intersection.insert(elem);
		}
	}
	std::cout << "标签集1大小为 " << largestComponent[0].size() << std::endl;
	std::cout << "标签集2大小为 " << largestComponent[1].size() << std::endl;
	std::cout << "最终交集大小为 " << intersection.size() << std::endl;

	if (intersection.find(querynode) == intersection.end())
	{
		return {};
	}

	tempgraph.filterEdgeLableAndRecover(labels);
	std::vector<std::unordered_set<int>> degreeVec = tempgraph.getOrderedNodes();
	std::unordered_map<int, int> nodeDegrees = tempgraph.getDegrees();
	std::unordered_set<int> community = result;

	int flag = false;

	int minDegree = 0;

	int nodeToDelete = -1;

	// 找到第一个待删除的节点
	for (int j = 0; j < k; j++)
	{
		if (!degreeVec[j].empty())
		{
			flag = true;
			minDegree = j;
			for (int node : degreeVec[j])
			{
				nodeToDelete = node;
				break;
			}
			break;
		}
	}
	while (flag)
	{
		// std::cout << community.size() << std::endl;
		// std::cout<<"the node to delete is "<<nodeToDelete<<std::endl;
		intersection.erase(nodeToDelete);
		degreeVec[nodeDegrees[nodeToDelete]].erase(nodeToDelete);
		for (int neighbor : tempgraph.getNeighbors(nodeToDelete))
		{
			if (intersection.find(neighbor) != intersection.end())
			{ // 仅更新仍在社区中的节点
				degreeVec[nodeDegrees[neighbor]].erase(neighbor);
				nodeDegrees[neighbor]--;
				degreeVec[nodeDegrees[neighbor]].insert(neighbor);
			}
		}
		flag = false;
		for (int j = 0; j < k; j++)
		{
			if (!degreeVec[j].empty())
			{
				flag = true;
				minDegree = j;
				for (int node : degreeVec[minDegree])
				{
					nodeToDelete = node;
					break;
				}
				break;
			}
		}
	} // 从查询节点开始检索最大连通分量

	return intersection;
}

std::unordered_set<int> Hierarchy::bl_index(std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &edgeLables, int querynode, int k, Graph &graph, std::unordered_set<std::string> &labels, LCS &lcs, bool isLCS)
{
	std::unordered_set<int> result;
	if (coreIndex[querynode] < k)
	{
		//std::cout << "init k-core unvaild"<<std::endl;
		return {};
	}
	if (isLCS&&!lcs.isExist(querynode, labels, k))
	{
		// result.insert(queryNodes);
		//std::cout << "init k-core unvaild"<<std::endl;
		return {};
	}
	Graph tempgraph(graph);
	if (!isLCS)
	{
		tempgraph.filterEdgeLableAndRecover(labels);
	}
	int targetNodeId = vertexToNode[querynode];

	// 简化层次遍历逻辑
	while (nodeCoreness[targetNodeId] >= k)
	{
		int parentNodeId = parents[targetNodeId];
		if (nodeCoreness[parentNodeId] < k || !isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
		{
			break;
		}
		targetNodeId = parentNodeId;
	}

	// 初始化两个标签对应的社区
	auto labelIter = labels.begin();
	std::string label1 = *labelIter++;
	std::string label2 = *labelIter;
	std::unordered_set<std::string> s1{label1};
	std::unordered_set<std::string> s2{label2};
	std::unordered_set<int> community1;
	std::unordered_set<int> community2;

	// BFS 遍历子树
	std::queue<int> queue;
	queue.push(targetNodeId);
	while (!queue.empty())
	{
		int current = queue.front();
		queue.pop();
		for (int vertex : nodes[current])
		{
			if (isLCS)
			{
				if (lcs.isExist(vertex, s1, k))
				{
					community1.insert(vertex);
				}
				if (lcs.isExist(vertex, s2, k))
				{
					community2.insert(vertex);
				}
			}
			else
			{
				community1.insert(vertex);
				community2.insert(vertex);
			}
		}
		for (int child : children[current])
		{
			queue.push(child);
		}
	}
	std::unordered_set<int> ans;
	// 处理第一个标签
	if (!isLCS)
	{
		std::unordered_set<std::string> s{label1};
		Graph temp(tempgraph);
		// for (int node : temp.getNodes())
		// {
		// 	if (community1.find(node) == community1.end())
		// 	{
		// 		temp.removeNode(node);
		// 	}
		// }

		temp.filterEdgeLableAndRecover(s1);
		ans = CommunitySearch::onlineWithEdge(temp, k, s);
		// 添加 BFS，获取包含查询节点的连通分量
		{
			std::unordered_set<int> visited;
			std::queue<int> q;
			q.push(querynode);
			visited.insert(querynode);

			while (!q.empty())
			{
				int current = q.front();
				q.pop();

				for (int neighbor : temp.getNeighbors(current))
				{
					if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end())
					{
						q.push(neighbor);
						visited.insert(neighbor);
					}
				}
			}
			ans = visited;
		}

		community1 = ans;
		//std::cout << "the community1 size is " << ans.size() << std::endl;
	}
	else
	{
		ans = community1;
		{
			std::unordered_set<int> visited;
			std::queue<int> q;
			q.push(querynode);
			visited.insert(querynode);

			while (!q.empty())
			{
				int current = q.front();
				q.pop();

				for (int neighbor : tempgraph.getNeighbors(current))
				{
					if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end() && edgeLables[current][neighbor].find(label1) != edgeLables[current][neighbor].end())
					{
						q.push(neighbor);
						visited.insert(neighbor);
					}
				}
			}
			ans = visited;
		}
		community1 = ans;
		//std::cout << "the community1 size is " << ans.size() << std::endl;
	}

	ans.clear();
	// 处理第2个标签
	if (!isLCS)
	{
		std::unordered_set<std::string> s{label2};
		// for (int node : temp.getNodes())
		// {
		// 	if (community2.find(node) == community2.end())
		// 	{
		// 		temp.removeNode(node);
		// 	}
		// }
		tempgraph.filterEdgeLableAndRecover(s2);
		ans = CommunitySearch::onlineWithEdge(tempgraph, k, s);
		// 添加 BFS，获取包含查询节点的连通分量

		{
			std::unordered_set<int> visited;
			std::queue<int> q;
			q.push(querynode);
			visited.insert(querynode);

			while (!q.empty())
			{
				int current = q.front();
				q.pop();

				for (int neighbor : tempgraph.getNeighbors(current))
				{
					if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end())
					{
						q.push(neighbor);
						visited.insert(neighbor);
					}
				}
			}
			// 更新 community2，只保留包含查询节点的连通分量
			ans = visited;
		}

		//std::cout << "the community2 size is " << ans.size() << std::endl;
		community2 = ans;
	}
	else
	{
		ans = community2;
		{
			std::unordered_set<int> visited;
			std::queue<int> q;
			q.push(querynode);
			visited.insert(querynode);

			while (!q.empty())
			{
				int current = q.front();
				q.pop();

				for (int neighbor : tempgraph.getNeighbors(current))
				{
					if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end() && edgeLables[current][neighbor].find(label2) != edgeLables[current][neighbor].end())
					{
						q.push(neighbor);
						visited.insert(neighbor);
					}
				}
			}
			ans = visited;
		}
		community2 = ans;
		//std::cout << "the community2 size is " << ans.size() << std::endl;
	}

	// 求交集
	//std::cout<<"the community1 size is "<<community1.size() << " the community2 size is " << community2.size() << std::endl;
	
	std::unordered_set<int> intersection;
	if(community1.size()==1&&community2.size()!=1){
		return community2;
	}else if (community1.size()!=1&&community2.size()==1){
		return community1;
	}

	for (const auto &elem : community1)
	{
		if (community2.find(elem) != community2.end())
		{
			intersection.insert(elem);
		}
	}

	return intersection;
}

std::unordered_set<int> Hierarchy::bl_mupt_index(std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &edgeLables, int querynode, int k, Graph &tempgraph, std::unordered_set<std::string> &labels, LCS &lcs, bool isLCS)
{
	if (!isLCS)
	{
		tempgraph.filterEdgeLableAndRecover(labels);
	}
	std::unordered_set<int> result;
	if (coreIndex[querynode] < k)
	{
		return {};
	}
	int targetNodeId = vertexToNode[querynode];

	// 简化层次遍历逻辑
	while (nodeCoreness[targetNodeId] >= k)
	{
		int parentNodeId = parents[targetNodeId];
		if (nodeCoreness[parentNodeId] < k || !isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
		{
			break;
		}
		targetNodeId = parentNodeId;
	}

	// 初始化标签对应的社区
	std::unordered_map<std::string, std::unordered_set<int>> community;

	// BFS 遍历子树
	std::queue<int> queue;
	queue.push(targetNodeId);
	while (!queue.empty())
	{
		int current = queue.front();
		queue.pop();
		for (int vertex : nodes[current])
		{
			if (isLCS)
			{
				for (std::string str : labels)
				{
					std::unordered_set<std::string> s{str};
					if (lcs.isExist(vertex, s, k))
					{
						community[str].insert(vertex);
					}
				}
			}
			else
			{
				for (std::string str : labels)
				{
					std::unordered_set<std::string> s{str};
					community[str].insert(vertex);
				}
			}
		}
		for (int child : children[current])
		{
			queue.push(child);
		}
	}
	std::unordered_set<int> ans;
	for (std::string label : labels)
	{
		if (!isLCS)
		{
			std::unordered_set<std::string> s{label};
			Graph temp(tempgraph);
			for (int node : temp.getNodes())
			{
				if (community[label].find(node) == community[label].end())
				{
					temp.removeNode(node);
				}
			}

			temp.filterEdgeLableAndRecover(s);
			ans = CommunitySearch::onlineWithEdge(temp, k, s);
			// 添加 BFS，获取包含查询节点的连通分量

			{
				std::unordered_set<int> visited;
				std::queue<int> q;
				q.push(querynode);
				visited.insert(querynode);

				while (!q.empty())
				{
					int current = q.front();
					q.pop();

					for (int neighbor : temp.getNeighbors(current))
					{
						if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end())
						{
							q.push(neighbor);
							visited.insert(neighbor);
						}
					}
				}
				ans = visited;
			}

			community[label] = ans;
			//std::cout << "the community1 size is " << ans.size() << std::endl;
		}
		else
		{
			ans = community[label];
			{
				std::unordered_set<int> visited;
				std::queue<int> q;
				q.push(querynode);
				visited.insert(querynode);

				while (!q.empty())
				{
					int current = q.front();
					q.pop();

					for (int neighbor : tempgraph.getNeighbors(current))
					{
						if (ans.find(neighbor) != ans.end() && visited.find(neighbor) == visited.end() && edgeLables[current][neighbor].find(label) != edgeLables[current][neighbor].end())
						{
							q.push(neighbor);
							visited.insert(neighbor);
						}
					}
				}
				ans = visited;
			}
			community[label] = ans;
			//std::cout << "the community1 size is " << ans.size() << std::endl;
		}
		ans.clear();
	}

	// 初始化交集集合，选择第一个集合作为初始集合
	auto it = community.begin();
	std::unordered_set<int> intersection = it->second;

	// 遍历 unordered_map 中其他的集合，并与当前交集集合进行取交集
	for (++it; it != community.end(); ++it)
	{
		std::unordered_set<int> tempIntersection;
		//std::cout<<"集合大小为 "<<it->second.size()<<" 当前标签为 "<<it->first<<std::endl;
		if (it->second.size() == 1)
			continue;
		if(intersection.size()==1&&it->second.size() != 1){
			intersection = it->second;
			continue;
		}
		for (const auto &elem : it->second)
		{
			if (intersection.find(elem) != intersection.end())
			{
				tempIntersection.insert(elem);
			}
		}
		intersection = std::move(tempIntersection);
	}
	//std::cout<<"交集大小为 "<<intersection.size()<<std::endl;

	return intersection;
}
// 两者是否有交集
bool Hierarchy::isSubsetOf(const std::unordered_set<std::string> &subset, const std::unordered_set<std::string> &set)
{
	for (const std::string &element : subset)
	{
		if (set.find(element) != set.end())
		{
			return true;
		}
	}
	for (const std::string &element : set)
	{
		if (subset.find(element) != subset.end())
		{
			return true;
		}
	}
	return false;
}

std::unordered_set<int> Hierarchy::queryWithLcs(Graph &graph, int queryNodes, int k, LCS &lcs, std::unordered_set<std::string> &labels)
{

	auto start1 = std::chrono::high_resolution_clock::now();

	std::unordered_set<int> result;
	std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &edgeLabels = graph.getEdgeLables();
	if (coreIndex[queryNodes] < k)
	{
		std::cout << "The query node is not in the k-core" << std::endl;
		return {};
	}
	if (!lcs.isExist(queryNodes, labels, k))
	{
		// result.insert(queryNodes);
		return {};
	}
	result.insert(queryNodes);
	int targetNodeId = vertexToNode[queryNodes];

	while (nodeCoreness[targetNodeId] >= k)
	{
		int parentNodeId = parents[targetNodeId];
		if (nodeCoreness[parentNodeId] <= k)
		{
			if (nodeCoreness[parentNodeId] < k)
			{
				break;
			}
			if (isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
			{
				targetNodeId = parentNodeId;
			}
			else
			{
				break;
			}
			// targetNodeId = parentNodeId;
		}
		else
		{
			if (isSubsetOf(edgeLabel[targetNodeId][parentNodeId], labels))
			{
				targetNodeId = parentNodeId;
			}
			else
			{
				break;
			}
		}
	}

	auto end1 = std::chrono::high_resolution_clock::now();
	auto start2 = std::chrono::high_resolution_clock::now();

	std::queue<int> queue;
	queue.push(targetNodeId);

	while (!queue.empty())
	{
		int current = queue.front();
		queue.pop();
		auto &currentNodes = nodes[current];
		for (int vertex : currentNodes)
		{
			if (lcs.isExist(vertex, labels, k))
			{
				result.insert(vertex);
			}
		}
		for (int child : children[current])
		{
			// if (isSubsetOf(edgeLabel[child][current], labels))
			// {
			// 	queue.push(child);
			// }else{
			// 	continue;
			// }
			queue.push(child);
		}
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	auto start3 = std::chrono::high_resolution_clock::now();

	std::unordered_set<int> resVertexs;
	std::unordered_set<int> visited;
	std::queue<int> tempQueue;
	tempQueue.push(queryNodes);
	resVertexs.insert(queryNodes);
	visited.insert(queryNodes);

	while (!tempQueue.empty())
	{
		int current = tempQueue.front();
		tempQueue.pop();
		for (int neighbor : graph.getNeighbors(current))
		{
			if (result.find(neighbor) != result.end() && visited.find(neighbor) == visited.end())
			{
				for (auto &label : edgeLabels[current][neighbor])
				{
					if (labels.find(label) != labels.end())
					{
						tempQueue.push(neighbor);
						resVertexs.insert(neighbor);
						visited.insert(neighbor);
					}
				}
			}
		}
	}
	auto end3 = std::chrono::high_resolution_clock::now();

	// std::cout << "Time for first loop: " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() << " ms.\n";
	// std::cout << "Time for second loop: " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << " ms.\n";
	// std::cout << "Time for third loop: " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count() << " ms.\n";

	return resVertexs;
}

int Hierarchy::printVerticesByCoreness()
{
	std::unordered_map<int, int> corenessToVertex;

	// 遍历 nodeCoreness 找到每层 coreness 对应的一个 vertexId
	for (const auto &pair : nodeCoreness)
	{
		int nodeId = pair.first;
		int coreness = pair.second;

		// 如果该 coreness 层还没有输出过 vertexId，就从 nodes 中取一个
		if (corenessToVertex.find(coreness) == corenessToVertex.end())
		{
			if (!nodes.at(nodeId).empty())
			{
				int vertexId = *nodes.at(nodeId).begin(); // 取当前 nodeId 下的一个 vertexId
				corenessToVertex[coreness] = vertexId;	  // 记录该 coreness 层的 vertexId
			}
		}
	}

	// 输出每层 coreness 对应的一个 vertexId
	for (const auto &pair : corenessToVertex)
	{
		// std::cout << "Coreness: " << pair.first << ", VertexId: " << pair.second << std::endl;
	}

	return kmax;
}

void Hierarchy::printData()
{
	std::cout << "Nodes:" << std::endl;
	for (const auto &pair : nodes)
	{
		std::cout << "Node " << pair.first << ": ";
		for (const auto &vertex : pair.second)
		{
			std::cout << vertex << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "Parents:" << std::endl;
	for (const auto &pair : parents)
	{
		std::cout << "Node " << pair.first << " has parent " << pair.second << std::endl;
	}

	std::cout << "Children:" << std::endl;
	for (const auto &pair : children)
	{
		std::cout << "Node " << pair.first << " has children: ";
		for (const auto &child : pair.second)
		{
			std::cout << child << " ";
		}
		std::cout << std::endl;
	}
}

std::unordered_map<int, std::unordered_set<int>> &Hierarchy::getCoreness()
{
	return coreness;
}

std::vector<int> Hierarchy::vaildCandicate(int i, int k)
{
	int K = 0.1 * i * k;
	std::vector<int> validNumbers;
	for (const auto &entry : coreness)
	{
		if (entry.first >= K)
		{
			for (const auto &num : entry.second)
			{
				validNumbers.push_back(num);
			}
		}
	}
	return validNumbers;
}