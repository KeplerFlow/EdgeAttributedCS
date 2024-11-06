#include "Graph.h"
#include <string>

int main(int argc, char *argv[])
{
    // Graph graph1("/home/asc23/lcy/graph/data/amazon/amazon-meta-graph-labeled-top8.txt");

    // Graph graph1("full.txt");

    // Graph graph1("/home/asc23/lcy/graph/Final/Dataset/movies_100k_network_full.txt");

    // Graph graph1("/home/asc23/lcy/graph/Final/Dataset/sem_network3_top8.txt");

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <graph file>" << std::endl;
        return 1;
    }

    std::string graphFile = argv[1];

    Graph graph1(graphFile);

    return 0;
}