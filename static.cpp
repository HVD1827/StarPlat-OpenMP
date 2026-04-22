#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <limits>
#include <fstream>
#include <chrono>

using namespace std;

struct Graph
{
    int n;
    vector<vector<pair<int, double>>> adj;

    Graph(int nodes = 0) : n(nodes), adj(nodes) {}

    void resize(int nodes)
    {
        n = nodes;
        adj.assign(nodes, {});
    }

    void add_edge(int u, int v, double w = 1.0)
    {
        if (u < 0 || v < 0)
            return;
        adj[u].push_back({v, w});
    }

    void add_undirected_edge(int u, int v, double w = 1.0)
    {
        add_edge(u, v, w);
        add_edge(v, u, w);
    }
};

// total edge weight / 2 for undirected
double total_edge_weight(const Graph &g)
{
    double total = 0.0;
    for (int i = 0; i < g.n; ++i)
        for (auto &e : g.adj[i])
            total += e.second;
    return total / 2.0;
}

double node_degree(const Graph &g, int node)
{
    double s = 0.0;
    for (auto &e : g.adj[node])
        s += e.second;
    return s;
}

// Phase 1: local moving
bool one_level(Graph &graph, vector<int> &nodeComm, double m)
{
    int n = graph.n;
    bool improvement = false;
    bool localMoved = true;

    vector<double> degrees(n);
    for (int i = 0; i < n; ++i)
        degrees[i] = node_degree(graph, i);

    unordered_map<int, double> commTot;
    for (int i = 0; i < n; ++i)
        commTot[nodeComm[i]] += degrees[i];

    while (localMoved)
    {
        localMoved = false;

        for (int node = 0; node < n; ++node)
        {
            int curComm = nodeComm[node];
            double k_i = degrees[node];

            unordered_map<int, double> k_i_in;
            for (auto &e : graph.adj[node])
            {
                int neigh = e.first;
                double w = e.second;
                int neighComm = nodeComm[neigh];
                if (neigh == node)
                    continue;
                k_i_in[neighComm] += w;
            }

            commTot[curComm] -= k_i;
            nodeComm[node] = -1;

            double bestGain = 0.0;
            int bestComm = curComm;

            for (auto &entry : k_i_in)
            {
                int targetComm = entry.first;
                double k_i_in_toC = entry.second;
                double totC = commTot[targetComm];
                double gain = (k_i_in_toC - (k_i * totC / (2.0 * m))) / (2.0 * m);
                if (gain > bestGain)
                {
                    bestGain = gain;
                    bestComm = targetComm;
                }
            }

            nodeComm[node] = bestComm;
            commTot[bestComm] += k_i;

            if (bestComm != curComm)
            {
                localMoved = true;
                improvement = true;
            }
        }
    }

    return improvement;
}

// Aggregate graph + mapping update
Graph aggregate_graph(const Graph &graph, const vector<int> &nodeComm,
                      vector<int> &newIndexOfComm)
{
    int n = graph.n;
    unordered_map<int, int> commMap;
    int idx = 0;
    for (int i = 0; i < n; ++i)
    {
        int c = nodeComm[i];
        if (commMap.find(c) == commMap.end())
            commMap[c] = idx++;
    }
    int newN = idx;
    newIndexOfComm.assign(newN, 0);

    map<pair<int, int>, double> edgeAcc;
    for (int u = 0; u < n; ++u)
    {
        int cu = commMap.at(nodeComm[u]);
        for (auto &e : graph.adj[u])
        {
            int v = e.first;
            int cv = commMap.at(nodeComm[v]);
            double w = e.second;
            pair<int, int> key = (cu <= cv) ? make_pair(cu, cv) : make_pair(cv, cu);
            edgeAcc[key] += w;
        }
    }

    Graph newG(newN);
    for (auto &entry : edgeAcc)
    {
        int a = entry.first.first;
        int b = entry.first.second;
        double w = entry.second;
        newG.add_undirected_edge(a, b, w / 2.0);
    }

    return newG;
}

// Louvain with global mapping
void louvain(Graph &graph, vector<int> &finalCommunities, int max_pass = 100)
{
    int n = graph.n;
    vector<int> nodeCommunity(n);
    for (int i = 0; i < n; ++i)
        nodeCommunity[i] = i;

    double m = total_edge_weight(graph);
    if (m <= 0)
        return;

    // ✅ Global mapping: maps original nodes to final communities
    vector<int> globalCommunity = nodeCommunity;

    bool improvement = true;
    int pass = 0;
    while (improvement && pass < max_pass)
    {
        pass++;
        bool moved = one_level(graph, nodeCommunity, m);
        if (!moved)
            break;

        vector<int> newIndexOfComm;
        Graph newG = aggregate_graph(graph, nodeCommunity, newIndexOfComm);

        // Build mapping from old community IDs → compact new ones
        unordered_map<int, int> commMap;
        int idx = 0;
        for (int c : nodeCommunity)
            if (commMap.find(c) == commMap.end())
                commMap[c] = idx++;

        // ✅ Update global mapping for all original nodes
        for (int i = 0; i < (int)globalCommunity.size(); ++i)
        {
            globalCommunity[i] = commMap[nodeCommunity[globalCommunity[i]]];
        }

        // prepare next iteration
        vector<int> newNodeComm(newG.n);
        for (int i = 0; i < newG.n; ++i)
            newNodeComm[i] = i;

        graph = std::move(newG);
        nodeCommunity = std::move(newNodeComm);
        m = total_edge_weight(graph);
    }

    finalCommunities = globalCommunity;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <graph_file>\n";
        return 1;
    }

    string graph_file = argv[1];

    ifstream fin(graph_file);
    if (!fin.is_open())
    {
        cerr << "Error opening graph file\n";
        return 1;
    }

    int n, m;
    fin >> n >> m;

    if (fin.fail())
    {
        cerr << "Error reading graph dimensions\n";
        return 1;
    }

    Graph g(n);

    for (int i = 0; i < m; i++)
    {
        int u, v;
        double w;
        fin >> u >> v >> w;

        if (fin.fail())
        {
            cerr << "Error reading edge " << i << "\n";
            return 1;
        }

        g.add_undirected_edge(u, v, w);
    }

    fin.close();

    auto start_static = chrono::high_resolution_clock::now();
    vector<int> communities;
    louvain(g, communities);
    auto end_static = chrono::high_resolution_clock::now();

    double static_time = chrono::duration<double>(end_static - start_static).count();
    cout << "\nInitial Louvain completed in " << static_time << " seconds.\n";

    // cout << "Initial communities:\n";
    // for (int i = 0; i < (int)communities.size(); ++i)
    //     cout << "Node " << i << " -> Community " << communities[i] << "\n";

    return 0;
}
