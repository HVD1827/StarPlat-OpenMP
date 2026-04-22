// louvain_fixed.cpp
#include <bits/stdc++.h>
#include <omp.h>
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
        // small critical to avoid concurrent push_back races if called concurrently
#pragma omp critical
        {
            add_edge(u, v, w);
            add_edge(v, u, w);
        }
    }

    void remove_undirected_edge(int u, int v)
    {
#pragma omp critical
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
    }
};

//
// ---------- Utility Functions ----------
//
double total_edge_weight(const Graph &g)
{
    double total = 0.0;
#pragma omp parallel for reduction(+ : total)
    for (int i = 0; i < g.n; ++i)
        for (auto &e : g.adj[i])
            total += e.second;
    return total / 2.0; // undirected
}

double node_degree(const Graph &g, int node)
{
    double s = 0.0;
    for (auto &e : g.adj[node])
        s += e.second;
    return s;
}

//
// ---------- One Level Optimization (SEQUENTIAL moves, parallel helpers) ----------
//
bool one_level(Graph &graph, vector<int> &nodeComm, double m)
{
    const int n = graph.n;
    bool improvement = false;
    bool localMoved = true;

    // parallel degrees
    vector<double> degrees(n, 0.0);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        degrees[i] = node_degree(graph, i);

    // compute community totals in parallel -> merged
    unordered_map<int, double> commTot;
#pragma omp parallel
    {
        unordered_map<int, double> local;
#pragma omp for nowait
        for (int i = 0; i < n; ++i)
            local[nodeComm[i]] += degrees[i];

#pragma omp critical
        for (auto &p : local)
            commTot[p.first] += p.second;
    }

    // Sequential node-move loop (safe & deterministic)
    while (localMoved)
    {
        localMoved = false;

        for (int node = 0; node < n; ++node)
        {
            int curComm = nodeComm[node];
            double k_i = degrees[node];

            // collect k_i_in (weights from node to each community)
            unordered_map<int, double> k_i_in;
            for (auto &e : graph.adj[node])
            {
                int neigh = e.first;
                if (neigh == node)
                    continue;
                k_i_in[nodeComm[neigh]] += e.second;
            }

            // temporarily remove node from its community totals
            commTot[curComm] -= k_i;
            nodeComm[node] = -1;

            // find best community to move into
            double bestGain = 0.0;
            int bestComm = curComm;
            for (auto &entry : k_i_in)
            {
                int target = entry.first;
                double k_i_in_toC = entry.second;
                double totC = commTot.count(target) ? commTot[target] : 0.0;
                double gain = (k_i_in_toC - (k_i * totC / (2.0 * m))) / (2.0 * m);
                if (gain > bestGain)
                {
                    bestGain = gain;
                    bestComm = target;
                }
            }

            // place node into best community (could be same as cur)
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
    commMap.reserve(n);
    int idx = 0;
    for (int i = 0; i < n; ++i)
    {
        int c = nodeComm[i];
        if (!commMap.count(c))
            commMap[c] = idx++;
    }

    // use local maps per thread then merge to avoid contention
    map<pair<int, int>, double> edgeAcc;
#pragma omp parallel
    {
        map<pair<int, int>, double> localAcc;
#pragma omp for nowait
        for (int u = 0; u < n; ++u)
        {
            int cu = commMap.at(nodeComm[u]);
            for (auto &e : graph.adj[u])
            {
                int v = e.first;
                int cv = commMap.at(nodeComm[v]);
                double w = e.second;
                pair<int, int> key = (cu <= cv) ? make_pair(cu, cv) : make_pair(cv, cu);
                localAcc[key] += w;
            }
        }
#pragma omp critical
        for (auto &p : localAcc)
            edgeAcc[p.first] += p.second;
    }

    Graph newG(idx);
    for (auto &entry : edgeAcc)
    {
        int a = entry.first.first;
        int b = entry.first.second;
        double w = entry.second;
        // each undirected edge has been counted twice in adjacency accumulation,
        // so divide by 2.0 when adding undirected edge weights
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
    iota(nodeComm.begin(), nodeComm.end(), 0);
    globalComm = nodeComm;

    double m = total_edge_weight(graph);
    if (m <= 0.0)
        return;

    bool improvement = true;
    int pass = 0;

    while (improvement && pass < max_pass)
    {
        pass++;
        bool moved = one_level(graph, nodeComm, m);
        if (!moved)
            break;

        // canonicalize community ids and write to globalComm for original nodes
        unordered_map<int, int> cmap;
        cmap.reserve(nodeComm.size());
        int idx = 0;
        for (int c : nodeComm)
            if (!cmap.count(c))
                cmap[c] = idx++;

        // If the original graph shrank in node count due to aggregation, we only
        // provide mapping for the current number of nodes (globalComm size will be updated below).
        // For simplicity, if we're at top-level we copy mapping back for nodes that still exist.
        if ((int)globalComm.size() == n)
        {
            for (int i = 0; i < n; ++i)
                globalComm[i] = cmap[nodeComm[i]];
        }
        else
        {
            // otherwise, just reset globalComm to nodeComm-like mapping of current size
            globalComm.resize(n);
            for (int i = 0; i < n; ++i)
                globalComm[i] = cmap[nodeComm[i]];
        }

        // aggregate graph to next level
        Graph newG = aggregate_graph(graph, nodeComm);

        // prepare nodeComm for next pass (one node per new community)
        nodeComm.assign(newG.n, 0);
        iota(nodeComm.begin(), nodeComm.end(), 0);

        graph = std::move(newG);
        m = total_edge_weight(graph);
    }
}

//
// ---------- Dynamic Updates ----------
//
void add_edge_dynamic(Graph &g, vector<int> &communities, int u, int v, double w = 1.0)
{
    // cout << "[Dynamic Update] Adding edge (" << u << ", " << v << ") w=" << w << "\n";
    g.add_undirected_edge(u, v, w);
    communities.clear();
    louvain(g, communities);
}

void remove_edge_dynamic(Graph &g, vector<int> &communities, int u, int v)
{
    // cout << "[Dynamic Update] Removing edge (" << u << ", " << v << ")\n";
    g.remove_undirected_edge(u, v);
    communities.clear();
    louvain(g, communities);
}

//
// ---------- Main ----------
//
int main(int argc, char **argv)
{
    // Optional: let user override thread count via OMP_NUM_THREADS or program arg

    if (argc < 5)
    {
        cout << "Usage: " << argv[0] << " <graph_file> <updates_file> <omp_threads> <num_updates>\n";
        return 1;
    }

    cout << atoi(argv[3]) << endl;
    int t = atoi(argv[3]);
    omp_set_num_threads(t);

    omp_set_num_threads(t);

#pragma omp parallel
    {
#pragma omp master
        {
            cout << "Number of threads: " << omp_get_num_threads() << endl;
        }
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
    if (!fin)
    {
        cerr << "Invalid graph header\n";
        return 1;
    }
    Graph g(n);

    for (int i = 0; i < m; ++i)
    {
        int u, v;
        double w;
        fin >> u >> v >> w;
        if (!fin)
        {
            cerr << "Invalid edge line\n";
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

    cout << "Initial communities size: " << communities.size() << "\n";
    for (size_t i = 0; i < communities.size() && i < 100; ++i)
        cout << "Node " << i << " -> Community " << communities[i] << "\n";

    ifstream upd(updates_file);
    if (!upd.is_open())
    {
        // it's okay to not have updates; just finish
        cout << "No updates file found; exiting after initial run.\n";
        return 0;
    }

    string line;
    int ln = 0;

    auto start_update = chrono::high_resolution_clock::now();

    int max_num_upd = atoi(argv[4]);
    int num_upd = 0;

    while (getline(upd, line))
    {
        num_upd++;
        if (num_upd > max_num_upd)
            break;

        ++ln;
        if (line.empty())
            continue;
        istringstream iss(line);
        char op;
        iss >> op;
        if (op == 'a')
        {
            int u, v;
            double w;
            if (!(iss >> u >> v >> w))
            {
                cerr << "Invalid add at line " << ln << ": " << line << "\n";
                continue;
            }
            add_edge_dynamic(g, communities, u, v, w);
        }
        else if (op == 'd')
        {
            int u, v;
            if (!(iss >> u >> v))
            {
                cerr << "Invalid del at line " << ln << ": " << line << "\n";
                continue;
            }
            remove_edge_dynamic(g, communities, u, v);
        }
        else
        {
            cerr << "Unknown op at line " << ln << ": " << line << "\n";
        }
    }
    upd.close();

    auto end_update = chrono::high_resolution_clock::now();

    double per_update_time = chrono::duration<double>(end_update - start_update).count();
    cout << "\nper update time " << per_update_time / (1.0 * ln) << " seconds.\n";

    cout << communities.size() << endl;
    // cout << "Final communities (first 100 shown if large):\n";
    // for (size_t i = 0; i < communities.size() && i < 100; ++i)
    //     cout << "Node " << i << " -> Community " << communities[i] << "\n";

    return 0;
}
