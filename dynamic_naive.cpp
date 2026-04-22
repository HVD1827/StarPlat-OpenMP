#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>
#include <chrono>

using namespace std;

//
// ---------- Graph Structure ----------
//
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
        if (u < 0 || v < 0 || u >= n || v >= n)
            return;
        adj[u].push_back({v, w});
    }

    void add_undirected_edge(int u, int v, double w = 1.0)
    {
        add_edge(u, v, w);
        add_edge(v, u, w);
    }

    void remove_undirected_edge(int u, int v)
    {
        auto remove_edge = [&](int a, int b)
        {
            auto &vec = adj[a];
            vec.erase(remove_if(vec.begin(), vec.end(),
                                [&](auto &p)
                                { return p.first == b; }),
                      vec.end());
        };
        remove_edge(u, v);
        remove_edge(v, u);
    }
};

//
// ---------- Utility Functions ----------
//
double total_edge_weight(const Graph &g)
{
    double total = 0.0;
    for (int i = 0; i < g.n; ++i)
        for (auto &e : g.adj[i])
            total += e.second;
    return total / 2.0; // because undirected
}

double node_degree(const Graph &g, int node)
{
    double s = 0.0;
    for (auto &e : g.adj[node])
        s += e.second;
    return s;
}

//
// ---------- One Level Optimization ----------
//
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

//
// ---------- Graph Aggregation ----------
//
Graph aggregate_graph(const Graph &graph, const vector<int> &nodeComm)
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

    map<pair<int, int>, double> edgeAcc;
    for (int u = 0; u < n; ++u)
    {
        int cu = commMap[nodeComm[u]];
        for (auto &e : graph.adj[u])
        {
            int v = e.first;
            int cv = commMap[nodeComm[v]];
            double w = e.second;
            pair<int, int> key = (cu <= cv) ? make_pair(cu, cv) : make_pair(cv, cu);
            edgeAcc[key] += w;
        }
    }

    Graph newG(idx);
    for (auto &entry : edgeAcc)
    {
        int a = entry.first.first;
        int b = entry.first.second;
        double w = entry.second;
        newG.add_undirected_edge(a, b, w / 2.0);
    }
    return newG;
}

//
// ---------- Louvain Algorithm ----------
//
void louvain(Graph graph, vector<int> &globalComm, int max_pass = 100)
{
    int n = graph.n;
    vector<int> nodeComm(n);
    for (int i = 0; i < n; ++i)
        nodeComm[i] = i;

    globalComm = nodeComm; // fresh initialization each run

    double m = total_edge_weight(graph);
    if (m <= 0)
        return;

    bool improvement = true;
    int pass = 0;

    while (improvement && pass < max_pass)
    {
        pass++;
        bool moved = one_level(graph, nodeComm, m);
        if (!moved)
            break;

        unordered_map<int, int> commMap;
        int idx = 0;
        for (int c : nodeComm)
            if (commMap.find(c) == commMap.end())
                commMap[c] = idx++;

        for (int i = 0; i < n; ++i)
            globalComm[i] = commMap[nodeComm[i]];

        graph = aggregate_graph(graph, nodeComm);
        vector<int> newNodeComm(graph.n);
        for (int i = 0; i < graph.n; ++i)
            newNodeComm[i] = i;
        nodeComm = move(newNodeComm);
        m = total_edge_weight(graph);
    }
}

//
// ---------- Dynamic Updates (Re-run Static) ----------
//
void add_edge_dynamic(Graph &g, vector<int> &communities, int u, int v, double w = 1.0)
{
    // cout << "\n[Dynamic Update] Adding edge (" << u << ", " << v << ")\n";
    g.add_undirected_edge(u, v, w);
    communities.clear(); // ensure fresh computation
    louvain(g, communities);
}

void remove_edge_dynamic(Graph &g, vector<int> &communities, int u, int v)
{
    // cout << "\n[Dynamic Update] Removing edge (" << u << ", " << v << ")\n";
    g.remove_undirected_edge(u, v);
    communities.clear(); // ensure fresh computation
    louvain(g, communities);
}

//
// ---------- Main ----------
//
int main(int argc, char *argv[])
{

    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " <graph_file> <updates_file> <num_updates>\n";
        return 1;
    }

    string graph_file = argv[1];
    string updates_file = argv[2];

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

    ifstream upd(updates_file);
    if (!upd.is_open())
    {
        cerr << "Error opening updates file\n";
        return 1;
    }

    string line;
    int lineNum = 0;

    auto start_update = chrono::high_resolution_clock::now();

    int max_num_upd = atoi(argv[3]);
    int num_upd = 0;

    while (getline(upd, line))
    {
        num_upd++;
        if (num_upd > max_num_upd)
            break;
        lineNum++;
        if (line.empty())
            continue; // Skip

        istringstream iss(line);
        char update;
        iss >> update;

        if (update == 'a')
        {
            int u, v;
            double w;
            if (iss >> u >> v >> w)
            {
                // cout << "Adding edge: " << u << " " << v << " " << w << "\n";
                add_edge_dynamic(g, communities, u, v, w);
            }
            else
            {
                cerr << "Error: Invalid ADD format on line " << lineNum << ": " << line << "\n";
            }
        }
        else if (update == 'd')
        {
            int u, v;
            if (iss >> u >> v)
            {
                // cout << "Removing edge: " << u << " " << v << "\n";
                remove_edge_dynamic(g, communities, u, v);
            }
            else
            {
                cerr << "Error: Invalid DELETE format on line " << lineNum << ": " << line << "\n";
            }
        }
        else
        {
            cerr << "Error: Unknown operation '" << update << "' on line " << lineNum << ": " << line << "\n";
        }
    }

    upd.close();

    auto end_update = chrono::high_resolution_clock::now();

    double per_update_time = chrono::duration<double>(end_update - start_update).count();
    cout << "\nper update time " << per_update_time / (1.0 * lineNum) << " seconds.\n";

    // cout << "\nFinal communities after updates:\n";
    // for (int i = 0; i < (int)communities.size(); ++i)
    //     cout << "Node " << i << " -> Community " << communities[i] << "\n";

    return 0;
}