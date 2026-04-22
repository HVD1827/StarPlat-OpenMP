#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <random>
#include <string>
#include <algorithm>

using namespace std;

// Function to read graph and generate updates
void generate_updates(const string& graph_file, const string& output_file, int num_updates, double add_ratio = 0.6) {
    ifstream fin(graph_file);
    if (!fin.is_open()) {
        cerr << "Error: Cannot open graph file " << graph_file << endl;
        return;
    }

    // Read graph dimensions
    int n, m;
    fin >> n >> m;

    // Read existing edges
    set<pair<int, int>> existing_edges;
    vector<pair<int, int>> edge_list;
    
    for (int i = 0; i < m; i++) {
        int u, v;
        double w;
        fin >> u >> v;
        if (fin.peek() != '\n' && fin.peek() != '\r') {
            fin >> w; // Read weight if present
        }
        
        // Store in canonical form (smaller node first)
        if (u > v) swap(u, v);
        existing_edges.insert({u, v});
        edge_list.push_back({u, v});
    }
    fin.close();

    cout << "Graph: " << n << " nodes, " << m << " edges" << endl;
    cout << "Generating " << num_updates << " updates (" << add_ratio * 100 << "% ADD, " << (1-add_ratio) * 100 << "% DELETE)" << endl;

    // Random number generation
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    uniform_int_distribution<> node_dist(0, n-1);

    vector<string> updates;
    int add_count = 0, delete_count = 0;
    int max_attempts = n * n; // Prevent infinite loops

    for (int i = 0; i < num_updates; i++) {
        if (dis(gen) < add_ratio) {
            // Generate ADD operation
            int attempts = 0;
            bool added = false;
            
            while (attempts < max_attempts && !added) {
                int u = node_dist(gen);
                int v = node_dist(gen);
                
                if (u != v) {
                    // Ensure canonical form
                    if (u > v) swap(u, v);
                    pair<int, int> edge = {u, v};
                    
                    if (existing_edges.find(edge) == existing_edges.end()) {
                        updates.push_back("a " + to_string(u) + " " + to_string(v) + " 1.0");
                        existing_edges.insert(edge);
                        edge_list.push_back(edge);
                        add_count++;
                        added = true;
                    }
                }
                attempts++;
            }
            
            if (!added) {
                cout << "Warning: Could not generate ADD operation after " << max_attempts << " attempts" << endl;
                // Skip this update or convert to DELETE
                i--;
            }
        } else {
            // Generate DELETE operation
            if (edge_list.empty()) {
                // No edges to delete, convert to ADD
                i--;
                continue;
            }
            
            uniform_int_distribution<> edge_dist(0, edge_list.size() - 1);
            int idx = edge_dist(gen);
            pair<int, int> edge = edge_list[idx];
            
            updates.push_back("d " + to_string(edge.first) + " " + to_string(edge.second));
            
            // Remove from both data structures
            existing_edges.erase(edge);
            edge_list.erase(edge_list.begin() + idx);
            delete_count++;
        }
    }

    // Write updates to file
    ofstream fout(output_file);
    if (!fout.is_open()) {
        cerr << "Error: Cannot open output file " << output_file << endl;
        return;
    }

    for (const string& update : updates) {
        fout << update << endl;
    }
    fout.close();

    cout << "Generated " << updates.size() << " updates:" << endl;
    cout << "  ADD operations: " << add_count << endl;
    cout << "  DELETE operations: " << delete_count << endl;
    cout << "  Remaining edges: " << existing_edges.size() << endl;
    cout << "Output written to: " << output_file << endl;
}

// Function to generate specific test cases
void generate_specific_updates(const string& graph_file, const string& output_file) {
    ifstream fin(graph_file);
    if (!fin.is_open()) {
        cerr << "Error: Cannot open graph file " << graph_file << endl;
        return;
    }

    int n, m;
    fin >> n >> m;
    fin.close();

    vector<string> updates = {
        "a 2 3 1.0",  // Connect two components
        "d 0 1",      // Remove bridge edge
        "a 1 4 1.0",  // Add cross-component edge
        "d 2 0",      // Remove another edge
        "a 0 3 1.0",  // Add another cross edge
        "d 1 2"       // Remove triangle edge
    };

    ofstream fout(output_file);
    for (const string& update : updates) {
        fout << update << endl;
    }
    fout.close();

    cout << "Generated " << updates.size() << " specific test updates" << endl;
    cout << "Output written to: " << output_file << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <graph_file> <output_file> [num_updates] [--specific]" << endl;
        cout << "  graph_file: Input graph file (n m, then edges u v [w])" << endl;
        cout << "  output_file: Output updates file" << endl;
        cout << "  num_updates: Number of updates to generate (default: 20)" << endl;
        cout << "  --specific: Generate specific test updates instead of random" << endl;
        return 1;
    }

    string graph_file = argv[1];
    string output_file = argv[2];
    int num_updates = 20;
    bool specific = false;

    // Parse command line arguments
    for (int i = 3; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--specific") {
            specific = true;
        } else {
            try {
                num_updates = stoi(arg);
            } catch (const exception& e) {
                cerr << "Error: Invalid number of updates '" << arg << "'" << endl;
                return 1;
            }
        }
    }

    if (specific) {
        generate_specific_updates(graph_file, output_file);
    } else {
        generate_updates(graph_file, output_file, num_updates, 0.6); // 60% adds
    }

    return 0;
}