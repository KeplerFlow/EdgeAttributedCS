#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>

#include "core.h"
#include "defs.h"
#include "gadget/gadget.h"
#include "glist/glist.h"
#include "traversal/traversal.h"
// 数据集的n为最大节点索引+1,m为(总行数-1)即0~m-1边
int main(int argc, char **argv)
{
  char path[100];     // path of input graph
  char method[100];   // method
  double ratio = 0.0; // ratio of edges to insert
  bool temp = false;  // temporal graph
  // initialize the options
  int option = -1;
  while (-1 != (option = getopt(argc, argv, "p:r:m:T:")))
  {
    switch (option)
    {
    case 'p':
      strcpy(path, optarg);
      break;
    case 'r':
      ASSERT(0.0 <= (ratio = atof(optarg)) && 1.0 >= ratio);
      break;
    }
  }
  // read the graph
  int n, m, m2;
  const auto read_func = gadget::ReadEdgesS;
  const auto edges = read_func(path, &n, &m);
  m2 = m - static_cast<int>(m * ratio);

  // initialize the core component
  core::GLIST *cm = new core::GLIST(n); // 创建原始对象
  core::GLISTData exportedData = cm->exportGLISTData();
  std::cout << "1";
  core::GLIST *cm2 = new core::GLIST(exportedData); // 使用拷贝构造函数创建副本
  std::cout << "12";
  // create the adjacent list representation
  std::vector<std::vector<int>> graph(n);
  for (int i = 0; i < m2; ++i)
  {
    int v1 = edges[i].first;
    int v2 = edges[i].second;
    graph[v1].push_back(v2);
    graph[v2].push_back(v1);
  }
  // compute the base core
  std::vector<int> core(n);

  cm->ComputeCore(graph, true, core);

  for (int i = m2; i < edges.size(); ++i)
  {
    cm->Insert(edges[i].first, edges[i].second, graph, core);
  }

  // verify the result
  {
    ERROR("check: insert", false);
    std::vector<int> tmp_core(n);

    cm->ComputeCore(graph, false, tmp_core);

    tmp_core[2] += 0;
    ASSERT_INFO(tmp_core == core, "wrong result after insert");
    cm->Check(graph, core);
    ERROR("check passed", false);
  }
  delete cm;
}
