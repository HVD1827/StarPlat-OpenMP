#include <bits/stdc++.h>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <input_graph> <output_graph>\n";
        return 1;
    }

    string input_graph = argv[1];
    string output_graph = argv[2];

    ifstream fin(input_graph);
    if (!fin.is_open())
    {
        cerr << "Error opening input file\n";
        return 1;
    }

    int n, m;
    fin >> n >> m;

    set<pair<int, int>> edges;
    for (int i = 0; i < m; ++i)
    {
        int u, v;
        fin >> u >> v;

        if (u == v)
            continue; // skip self-loops

        if (u > v)
            swap(u, v); // store in sorted order

        edges.insert({u, v}); // duplicates auto-removed
    }

    fin.close();

    ofstream fout(output_graph);
    if (!fout.is_open())
    {
        cerr << "Error opening output file\n";
        return 1;
    }

    fout << n << " " << edges.size() << "\n";
    for (auto &e : edges)
    {
        fout << e.first << " " << e.second << " 1.0\n";
    }

    fout.close();

    cout << "Filtered graph written to " << output_graph << "\n";
    cout << "Unique undirected edges: " << edges.size() << "\n";

    return 0;
}
