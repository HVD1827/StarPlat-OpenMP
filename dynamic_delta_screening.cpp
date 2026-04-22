// louvain_delta_fixed_full_fixed.cpp
#include <bits/stdc++.h>
#include <chrono>

using namespace std;

const double EPS_GAIN = 1e-12;

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
        auto rem = [&](int a, int b)
        {
            if (a < 0 || a >= n)
                return;
            auto &vec = adj[a];
            vec.erase(remove_if(vec.begin(), vec.end(), [&](auto &p)
                                { return p.first == b; }),
                      vec.end());
        };
        rem(u, v);
        rem(v, u);
    }
};

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

// compute modularity Q (direct formula) for verification
double compute_modularity(const Graph &g, const vector<int> &comm, const vector<double> &degrees)
{
    double m = total_edge_weight(g);
    if (m <= 0.0)
        return 0.0;
    double sum = 0.0;
    for (int i = 0; i < g.n; ++i)
    {
        for (auto &e : g.adj[i])
        {
            int j = e.first;
            if (comm[i] == comm[j])
            {
                sum += e.second - (degrees[i] * degrees[j]) / (2.0 * m);
            }
        }
    }
    return sum / (2.0 * m); // adjacency double-counts undirected edges
}

// --- helper: safe map increment ---
static inline void map_add(unordered_map<int, double> &m, int key, double val)
{
    auto it = m.find(key);
    if (it == m.end())
        m.emplace(key, val);
    else
        it->second += val;
}

// --- Static Louvain (single level, no aggregation) ---
// This implementation:
//  - considers neighbor communities (k_i_in) and also the singleton candidate (community id = node)
//  - uses unordered_map<int,double> for comm totals so new community ids can be created (we use node index as a possible new id)
void louvain_static(Graph &graph, vector<int> &nodeCommunity, int max_pass = 100)
{
    int n = graph.n;
    nodeCommunity.resize(n);
    for (int i = 0; i < n; ++i)
        nodeCommunity[i] = i; // each node its own community initially

    double m = total_edge_weight(graph);
    if (m <= 0.0)
        return;

    vector<double> degrees(n);
    for (int i = 0; i < n; ++i)
        degrees[i] = node_degree(graph, i);

    bool any_improvement = true;
    int pass = 0;
    while (any_improvement && pass < max_pass)
    {
        ++pass;
        any_improvement = false;

        // build commTot map (sum of degrees per community label)
        unordered_map<int, double> commTot;
        for (int i = 0; i < n; ++i)
            map_add(commTot, nodeCommunity[i], degrees[i]);

        bool localMoved = true;
        while (localMoved)
        {
            localMoved = false;
            for (int node = 0; node < n; ++node)
            {
                int curComm = nodeCommunity[node];
                double k_i = degrees[node];

                // compute k_i_in per community among neighbors
                unordered_map<int, double> k_i_in;
                for (auto &e : graph.adj[node])
                {
                    int neigh = e.first;
                    if (neigh == node)
                        continue;
                    int neighComm = nodeCommunity[neigh];
                    map_add(k_i_in, neighComm, e.second);
                }

                // remove node's contribution from its community totals temporarily
                commTot[curComm] -= k_i;

                // Evaluate best candidate: all neighbor communities + singleton (community id = node)
                double bestGain = 0.0;
                int bestComm = curComm;

                // evaluate neighbor communities
                for (auto &p : k_i_in)
                {
                    int target = p.first;
                    double kin = p.second;
                    double totC = 0.0;
                    auto it = commTot.find(target);
                    if (it != commTot.end())
                        totC = it->second;
                    // standard simplified Louvain gain:
                    double gain = (kin - (k_i * totC / (2.0 * m))) / (2.0 * m);
                    if (gain > bestGain + EPS_GAIN)
                    {
                        bestGain = gain;
                        bestComm = target;
                    }
                }

                // evaluate singleton (new community where totC = 0 and kin = 0 unless neighbor links to itself)
                // We consider candidate community id = node (unique per node). totC_new = commTot[node] (likely 0 unless reused), but safe to read.
                {
                    double kin = 0.0; // by construction a new empty community has no internal links from node (unless self-loop)
                    // however if some nodes already use community label == node, totC may be > 0; that's acceptable: we treat label=node as a candidate
                    double totC = 0.0;
                    auto it = commTot.find(node);
                    if (it != commTot.end())
                        totC = it->second;
                    double gain = (kin - (k_i * totC / (2.0 * m))) / (2.0 * m);
                    if (gain > bestGain + EPS_GAIN)
                    {
                        bestGain = gain;
                        bestComm = node;
                    }
                }

                // move if better
                if (bestComm != curComm)
                {
                    nodeCommunity[node] = bestComm;
                    localMoved = true;
                    any_improvement = true;
                }

                // add node back to its chosen community totals
                commTot[nodeCommunity[node]] += k_i;
            }
        }
    }
}

// --- One-level delta-screening active-set variant (uses a state of node contributions + active flags) ---
// This version will consider neighbor communities and the singleton candidate (node index)
struct DeltaScreeningState
{
    vector<double> nodeModularityContrib;
    vector<char> activeNodes;
    double threshold;
    DeltaScreeningState(int n = 0, double thr = 1e-6) : nodeModularityContrib(n, 0.0), activeNodes(n, 0), threshold(thr) {}
};

// calculate node modularity contribution (same simplified measure used in screening)
double calculate_node_modularity_contribution(const Graph &graph, int node,
                                              const vector<int> &communities,
                                              double m, const vector<double> &degrees,
                                              const unordered_map<int, double> *commTotPtr = nullptr)
{
    if (m <= 0.0)
        return 0.0;
    int nodeComm = communities[node];
    double k_i = degrees[node];
    double k_i_in = 0.0;
    for (auto &e : graph.adj[node])
    {
        int neigh = e.first;
        if (neigh == node)
            continue;
        if (communities[neigh] == nodeComm)
            k_i_in += e.second;
    }

    double totC = 0.0;
    if (commTotPtr)
    {
        auto it = commTotPtr->find(nodeComm);
        if (it != commTotPtr->end())
            totC = it->second;
        else
            totC = 0.0;
    }
    else
    {
        for (int i = 0; i < graph.n; ++i)
            if (communities[i] == nodeComm)
                totC += degrees[i];
    }

    double gain_num = k_i_in - (k_i * totC / (2.0 * m));
    return gain_num / (2.0 * m);
}

void mark_affected_nodes_delta(const Graph &graph, DeltaScreeningState &state,
                               int u, int v, const vector<int> &communities,
                               double m, const vector<double> &degrees)
{
    if ((int)state.nodeModularityContrib.size() != graph.n)
    {
        state = DeltaScreeningState(graph.n, state.threshold);
        unordered_map<int, double> commTot;
        for (int i = 0; i < graph.n; ++i)
            map_add(commTot, communities[i], degrees[i]);
        for (int i = 0; i < graph.n; ++i)
            state.nodeModularityContrib[i] = calculate_node_modularity_contribution(graph, i, communities, m, degrees, &commTot);
    }

    if (u >= 0 && u < graph.n)
        state.activeNodes[u] = 1;
    if (v >= 0 && v < graph.n)
        state.activeNodes[v] = 1;

    unordered_map<int, double> commTot;
    for (int i = 0; i < graph.n; ++i)
        map_add(commTot, communities[i], degrees[i]);

    if (u >= 0 && u < graph.n)
        state.nodeModularityContrib[u] = calculate_node_modularity_contribution(graph, u, communities, m, degrees, &commTot);
    if (v >= 0 && v < graph.n)
        state.nodeModularityContrib[v] = calculate_node_modularity_contribution(graph, v, communities, m, degrees, &commTot);

    auto check_neighbors = [&](int node)
    {
        if (node < 0 || node >= graph.n)
            return;
        for (auto &e : graph.adj[node])
        {
            int neigh = e.first;
            double oldMod = state.nodeModularityContrib[neigh];
            double newMod = calculate_node_modularity_contribution(graph, neigh, communities, m, degrees, &commTot);
            if (fabs(newMod - oldMod) > state.threshold)
                state.activeNodes[neigh] = 1;
        }
    };

    check_neighbors(u);
    check_neighbors(v);
}

bool one_level_delta_screening(Graph &graph, vector<int> &nodeComm, double m, DeltaScreeningState &state)
{
    int n = graph.n;
    if (n == 0 || m <= 0.0)
        return false;
    bool improvement = false;

    vector<double> degrees(n);
    for (int i = 0; i < n; ++i)
        degrees[i] = node_degree(graph, i);
    unordered_map<int, double> commTot;
    for (int i = 0; i < n; ++i)
        map_add(commTot, nodeComm[i], degrees[i]);

    unordered_set<int> activeSet;
    for (int i = 0; i < n; ++i)
        if (state.activeNodes[i])
            activeSet.insert(i);

    while (!activeSet.empty())
    {
        unordered_set<int> nextActive;
        for (int node : activeSet)
        {
            if (node < 0 || node >= n)
                continue;
            int curComm = nodeComm[node];
            double k_i = degrees[node];

            unordered_map<int, double> k_i_in;
            for (auto &e : graph.adj[node])
            {
                int neigh = e.first;
                if (neigh == node)
                    continue;
                int neighComm = nodeComm[neigh];
                map_add(k_i_in, neighComm, e.second);
            }

            // remove node temporarily from curComm totals
            commTot[curComm] -= k_i;

            double bestGain = 0.0;
            int bestComm = curComm;

            // consider neighbor communities
            for (auto &p : k_i_in)
            {
                int target = p.first;
                double kin = p.second;
                double totC = 0.0;
                auto it = commTot.find(target);
                if (it != commTot.end())
                    totC = it->second;
                double gain = (kin - (k_i * totC / (2.0 * m))) / (2.0 * m);
                if (gain > bestGain + EPS_GAIN)
                {
                    bestGain = gain;
                    bestComm = target;
                }
            }

            // consider singleton candidate (use label == node as candidate)
            {
                double kin = 0.0;
                double totC = 0.0;
                auto it = commTot.find(node);
                if (it != commTot.end())
                    totC = it->second;
                double gain = (kin - (k_i * totC / (2.0 * m))) / (2.0 * m);
                if (gain > bestGain + EPS_GAIN)
                {
                    bestGain = gain;
                    bestComm = node;
                }
            }

            // assign node to bestComm
            if (bestComm != curComm)
            {
                nodeComm[node] = bestComm;
                improvement = true;

                // update stored contribution for node
                state.nodeModularityContrib[node] = calculate_node_modularity_contribution(graph, node, nodeComm, m, degrees, &commTot);

                // neighbors might change significantly
                for (auto &e : graph.adj[node])
                {
                    int neigh = e.first;
                    double oldMod = state.nodeModularityContrib[neigh];
                    double newMod = calculate_node_modularity_contribution(graph, neigh, nodeComm, m, degrees, &commTot);
                    if (fabs(newMod - oldMod) > state.threshold)
                    {
                        nextActive.insert(neigh);
                        state.nodeModularityContrib[neigh] = newMod;
                    }
                }
            }

            // add node's degree back to (possibly new) community total
            commTot[nodeComm[node]] += k_i;
        }

        activeSet = move(nextActive);
    }

    // clear active flags
    fill(state.activeNodes.begin(), state.activeNodes.end(), 0);
    return improvement;
}

// --- dynamic add/remove wrappers that recompute baseline contributions safely ---
void add_edge_delta_screening(Graph &g, vector<int> &communities, DeltaScreeningState &state, int u, int v, double w = 1.0)
{
    // cout << "[Delta-Screening] Adding edge (" << u << "," << v << ")\n";
    if (u < 0 || v < 0 || u >= g.n || v >= g.n)
    {
        cerr << "add_edge indices out of range\n";
        return;
    }
    g.add_undirected_edge(u, v, w);
    double m_new = total_edge_weight(g);

    if ((int)communities.size() != g.n)
        louvain_static(g, communities);

    vector<double> degrees(g.n);
    for (int i = 0; i < g.n; ++i)
        degrees[i] = node_degree(g, i);

    if ((int)state.nodeModularityContrib.size() != g.n)
        state = DeltaScreeningState(g.n, state.threshold);

    unordered_map<int, double> commTot;
    for (int i = 0; i < g.n; ++i)
        map_add(commTot, communities[i], degrees[i]);
    for (int i = 0; i < g.n; ++i)
        state.nodeModularityContrib[i] = calculate_node_modularity_contribution(g, i, communities, m_new, degrees, &commTot);

    mark_affected_nodes_delta(g, state, u, v, communities, m_new, degrees);
    bool changed = one_level_delta_screening(g, communities, m_new, state);
    // if (changed)
    //     cout << "[Delta-Screening] Communities changed after addition.\n";
    // else
    //     cout << "[Delta-Screening] No community change after addition.\n";
}

void remove_edge_delta_screening(Graph &g, vector<int> &communities, DeltaScreeningState &state, int u, int v)
{
    // cout << "[Delta-Screening] Removing edge (" << u << "," << v << ")\n";
    if (u < 0 || v < 0 || u >= g.n || v >= g.n)
    {
        cerr << "remove_edge indices out of range\n";
        return;
    }
    g.remove_undirected_edge(u, v);
    double m_new = total_edge_weight(g);

    if ((int)communities.size() != g.n)
        louvain_static(g, communities);

    vector<double> degrees(g.n);
    for (int i = 0; i < g.n; ++i)
        degrees[i] = node_degree(g, i);

    if ((int)state.nodeModularityContrib.size() != g.n)
        state = DeltaScreeningState(g.n, state.threshold);

    unordered_map<int, double> commTot;
    for (int i = 0; i < g.n; ++i)
        map_add(commTot, communities[i], degrees[i]);
    for (int i = 0; i < g.n; ++i)
        state.nodeModularityContrib[i] = calculate_node_modularity_contribution(g, i, communities, m_new, degrees, &commTot);

    mark_affected_nodes_delta(g, state, u, v, communities, m_new, degrees);
    bool changed = one_level_delta_screening(g, communities, m_new, state);
    // if (changed)
    //     cout << "[Delta-Screening] Communities changed after removal.\n";
    // else
    //     cout << "[Delta-Screening] No community change after removal.\n";
}

// --- main test driver ---
int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " <graph_file> <updates_file> <num_updates>\n";
        return 1;
    }

    string graph_file = argv[1];
    string updates_file = argv[2];

    // ===== Load Graph =====
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

    // ===== Initial Louvain =====
    vector<int> communities;
    DeltaScreeningState state(n, 1e-6);

    auto start_static = chrono::high_resolution_clock::now();
    louvain_static(g, communities);
    auto end_static = chrono::high_resolution_clock::now();

    double static_time = chrono::duration<double>(end_static - start_static).count();
    cout << "\nInitial Louvain completed in " << static_time << " seconds.\n";

    // ===== Initialize Modularity Structures =====
    double m_total = total_edge_weight(g);
    if (m_total <= 0.0)
    {
        cout << "Empty or zero-weight graph\n";
        return 0;
    }

    vector<double> degrees(g.n);
    for (int i = 0; i < g.n; ++i)
        degrees[i] = node_degree(g, i);

    state = DeltaScreeningState(g.n, state.threshold);
    unordered_map<int, double> commTot;
    for (int i = 0; i < g.n; ++i)
        map_add(commTot, communities[i], degrees[i]);
    for (int i = 0; i < g.n; ++i)
        state.nodeModularityContrib[i] = calculate_node_modularity_contribution(g, i, communities, m_total, degrees, &commTot);

    // cout << "Initial communities:\n";
    // for (int i = 0; i < (int)communities.size(); ++i)
    //     cout << "Node " << i << " -> Community " << communities[i] << "\n";
    // cout << "Initial modularity Q = " << compute_modularity(g, communities, degrees) << "\n\n";

    // ===== Process Updates =====
    ifstream upd(updates_file);
    if (!upd.is_open())
    {
        cerr << "Error opening updates file\n";
        return 1;
    }

    cout << "Processing updates...\n";
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
            continue;

        istringstream iss(line);
        char type;
        iss >> type;

        if (type == 'A' || type == 'a')
        {
            int u, v;
            double w;
            if (iss >> u >> v >> w)
            {
                // cout << "Adding edge: " << u << " " << v << " " << w << "\n";
                add_edge_delta_screening(g, communities, state, u, v, w);
            }
            else
            {
                cerr << "Error: Invalid ADD format on line " << lineNum << ": " << line << "\n";
            }
        }
        else if (type == 'D' || type == 'd')
        {
            int u, v;
            if (iss >> u >> v)
            {
                // cout << "Removing edge: " << u << " " << v << "\n";
                remove_edge_delta_screening(g, communities, state, u, v);
            }
            else
            {
                cerr << "Error: Invalid DELETE format on line " << lineNum << ": " << line << "\n";
            }
        }
        else
        {
            cerr << "Error: Unknown operation '" << type << "' on line " << lineNum << ": " << line << "\n";
        }
    }

    upd.close();

    auto end_update = chrono::high_resolution_clock::now();
    double total_update_time = chrono::duration<double>(end_update - start_update).count();

    cout << "\nTotal update time: " << total_update_time << " seconds.\n";
    if (lineNum > 0)
        cout << "Average per-update time: " << (total_update_time / lineNum) << " seconds.\n";

    // ===== Final Output =====
    // cout << "\nFinal communities after updates:\n";
    // for (int i = 0; i < (int)communities.size(); ++i)
    //     cout << "Node " << i << " -> Community " << communities[i] << "\n";

    // for (int i = 0; i < g.n; ++i)
    //     degrees[i] = node_degree(g, i);
    // cout << "Final modularity Q = " << compute_modularity(g, communities, degrees) << "\n";

    return 0;
}
