#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstring>
#include <climits>
#include "/home/harsh/Desktop/Sem7/UGRC-I/starplat/graphcode/generated_omp/dynamicBatchTCV2_dyn.cc"

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <graphFile> <updateFile>\n";
        return 1;
    }

    char *graphFile = argv[1];
    char *updateFile = argv[2];

    graph g(graphFile);
    g.parseGraph();
    std::cout << "Initial triangles: " << staticTC(g) << std::endl;

    std::vector<update> updates = g.parseUpdates(updateFile);

    int batchSize = 10;
    DynTC(g, updates, batchSize);

    return 0;
}
