#ifndef HIERARCHY_H
#define HIERARCHY_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <climits>
#include <chrono>
#include <cassert>
#include <random>
#include <ctime>

#include "../Graph.h"
#include "../Algorithm/CoreGroup.h"
#include "../Algorithm/Online.h"
#include "LCS.h"

class Hierarchy
{
private:
	// key is the nodeId , value is the vertexs;
	std::unordered_map<int, std::unordered_set<int>> nodes;
	// key is the nodeId , value is the coreness;
	std::unordered_map<int, int> nodeCoreness;
	// key is the vertexId , value is the nodeId;
	std::unordered_map<int, int> vertexToNode;
	// key is the nodeId , value is the parent nodeId;
	std::unordered_map<int, int> parents;
	// key is the nodeId , value is the children nodeId;
	std::unordered_map<int, std::unordered_set<int>> children;

private:
	// 建立索引保存父索引的边标签,包括边标签,顶点对,树形索引节点对
	std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> edgeLabel;
	// key is the vertex , value is the coreness;
	std::unordered_map<int, int> coreIndex;
	// key is the coreness, vlaue is the vertex set;
	std::unordered_map<int, std::unordered_set<int>> coreness;
	// record the nodes' scale
	int nodeSum;
	int kmax;
	int vertexSum;

public:
	Hierarchy(Graph &graph);
	int coreDecomposition(Graph &graph);
	std::unordered_map<int, std::unordered_set<int>> BFS(Graph &graph, int startVertex, std::unordered_map<int, std::unordered_set<int>> &vertexs, int minCoreness, int nodeId, int &curNodeId);
	void DFS(Graph &graph, std::unordered_map<int, std::unordered_set<int>> &vertexs, int minCoreness, int nodeId);
	void build(Graph &graph);
	std::unordered_map<int, std::unordered_set<int>> getNodes() const { return nodes; }
	std::unordered_set<int> query(int queryNodes, int k);
	std::unordered_set<int> bl(int querynode, int k, Graph &graph, std::unordered_set<std::string> &labels);
	std::unordered_set<int> bl_mupt_index(std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &edgeLables, int querynode, int k, Graph &graph, std::unordered_set<std::string> &labels, LCS &lcs, bool isLCS);
	std::unordered_set<int> bl_index(std::unordered_map<int, std::unordered_map<int, std::unordered_set<std::string>>> &edgeLables, int querynode, int k, Graph &graph, std::unordered_set<std::string> &labels, LCS &lcs, bool isLCS);
	std::unordered_set<int> queryWithLcs(Graph &graph, int queryNodes, int k, LCS &lcs, std::unordered_set<std::string> &labels);

	void printData();
	bool isSubsetOf(const std::unordered_set<std::string> &subset, const std::unordered_set<std::string> &set);
	int printVerticesByCoreness();
	std::unordered_map<int, std::unordered_set<int>> &getCoreness();
	std::vector<int> vaildCandicate(int i, int k);
};

#endif