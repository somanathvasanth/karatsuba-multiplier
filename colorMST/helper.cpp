#include <bits/stdc++.h>
using namespace std;

typedef pair<int, int> iPair;
vector<vector<pair<int, int>>> adjMatrix;

struct Graph {
    int V, E;
    vector<pair<pair<int, int>, iPair>> edges;

    Graph(int V, int E) : V(V), E(E) {}

    void addEdge(int u, int v, int w, int red) {
        edges.push_back({{w, red}, {u, v}});
    }
};

struct DisjointSets {
    int *parent, *rnk;
    int n;

    DisjointSets(int n) : n(n) {
        parent = new int[n + 1];
        rnk = new int[n + 1];
        for (int i = 0; i <= n; i++) {
            rnk[i] = 0;
            parent[i] = i;
        }
    }

    int find(int u) {
        if (u != parent[u]) parent[u] = find(parent[u]);
        return parent[u];
    }

    void merge(int x, int y) {
        x = find(x), y = find(y);
        if (rnk[x] > rnk[y]) parent[y] = x;
        else parent[x] = y;
        if (rnk[x] == rnk[y]) rnk[y]++;
    }

    ~DisjointSets() {
        delete[] parent;
        delete[] rnk;
    }
};

pair<int, vector<pair<pair<int, int>, iPair>>> kruskalMST(Graph& g) {
    int mst_wt = 0;
    vector<pair<pair<int, int>, iPair>> mst_edges;
    sort(g.edges.begin(), g.edges.end());
    DisjointSets ds(g.V);

    for (auto& edge : g.edges) {
        int u = edge.second.first;
        int v = edge.second.second;
        int set_u = ds.find(u);
        int set_v = ds.find(v);
        if (set_u != set_v) {
            mst_edges.push_back(edge);
            mst_wt += edge.first.first;
            ds.merge(set_u, set_v);
        }
    }
    return {mst_wt, mst_edges};
}

pair<vector<int>, vector<int>> dfs(const vector<vector<pair<int, int>>>& adjMatrix, int startVertex) {
    int V = adjMatrix.size();
    vector<int> parent(V, -1), depth(V, -1);
    vector<bool> visited(V, false);
    stack<int> s;
    s.push(startVertex);
    visited[startVertex] = true;
    parent[startVertex] = startVertex;
    depth[startVertex] = 0;

    while (!s.empty()) {
        int u = s.top();
        s.pop();
        for (int v = 0; v < V; ++v) {
            if (adjMatrix[u][v].first != -1 && !visited[v]) {
                visited[v] = true;
                parent[v] = u;
                depth[v] = depth[u] + 1;
                s.push(v);
            }
        }
    }
    return {parent, depth};
}

vector<pair<pair<int, int>, iPair>> findCycleEdges(const vector<int>& parent, const vector<int>& depth, const pair<pair<int, int>, iPair>& new_edge) {
    int v1 = new_edge.second.first, v2 = new_edge.second.second;
    vector<pair<pair<int, int>, iPair>> cycle_edges;
    while (depth[v1] > depth[v2]) {
        cycle_edges.push_back({adjMatrix[v1][parent[v1]], {v1, parent[v1]}});
        v1 = parent[v1];
    }
    while (depth[v2] > depth[v1]) {
        cycle_edges.push_back({adjMatrix[v2][parent[v2]], {v2, parent[v2]}});
        v2 = parent[v2];
    }
    while (v1 != v2) {
        cycle_edges.push_back({adjMatrix[v1][parent[v1]], {v1, parent[v1]}});
        cycle_edges.push_back({adjMatrix[v2][parent[v2]], {v2, parent[v2]}});
        v1 = parent[v1];
        v2 = parent[v2];
    }
    return cycle_edges;
}

pair<pair<int, int>, iPair> findMaxWeightRedEdge(const vector<pair<pair<int, int>, iPair>>& cycle_edges) {
    pair<pair<int, int>, iPair> max_red_edge = {{-1, 1}, {-1, -1}};
    for (const auto& edge : cycle_edges) {
        if (edge.first.second == 1 && edge.first.first > max_red_edge.first.first) {
            max_red_edge = edge;
        }
    }
    return max_red_edge;
}

int main() {
    int V, E, threshold;
    cin >> V >> E >> threshold;
    Graph g(V, E);

    for (int i = 0; i < E; i++) {
        int u, v, w, r;
        cin >> u >> v >> w >> r;
        g.addEdge(u, v, w, r);
    }

    auto [mst_wt, mst_edges] = kruskalMST(g);
    adjMatrix.resize(V, vector<pair<int, int>>(V, {-1, -1}));
    for (const auto& edge : mst_edges) {
        int u = edge.second.first, v = edge.second.second;
        adjMatrix[u][v] = adjMatrix[v][u] = {edge.first.first, edge.first.second};
    }

    vector<pair<pair<int, int>, iPair>> blue_edges;
    for (const auto& edge : g.edges) {
        if (edge.first.second == 0) {
            bool in_mst = false;
            for (const auto& mst_edge : mst_edges) {
                if (edge == mst_edge) {
                    in_mst = true;
                    break;
                }
            }
            if (!in_mst) blue_edges.push_back(edge);
        }
    }
    sort(blue_edges.begin(), blue_edges.end());

    auto [parent, depth] = dfs(adjMatrix, 0);
    int red_edge_count = count_if(mst_edges.begin(), mst_edges.end(), [](auto e) { return e.first.second == 1; });

    while (!blue_edges.empty()) {
        // Corrected type: {weight_diff, {blue_edge, red_edge}}
        vector<pair<int, pair<pair<pair<int, int>, iPair>, pair<pair<int, int>, iPair>>>> candidates;
        for (auto& blue_edge : blue_edges) {
            int u = blue_edge.second.first, v = blue_edge.second.second;
            if (adjMatrix[u][v].first != -1) continue;

            auto cycle_edges = findCycleEdges(parent, depth, blue_edge);
            auto max_red_edge = findMaxWeightRedEdge(cycle_edges);
            if (max_red_edge.second.first != -1) {
                int new_wt = mst_wt + blue_edge.first.first - max_red_edge.first.first;
                if (new_wt <= threshold) {
                    candidates.push_back({blue_edge.first.first - max_red_edge.first.first, {blue_edge, max_red_edge}});
                }
            }
        }

        if (candidates.empty()) break;

        auto best = *min_element(candidates.begin(), candidates.end()); // Minimize weight increase
        auto blue_edge = best.second.first;  // Correctly extract blue_edge
        auto red_edge = best.second.second;  // Correctly extract red_edge

        // Update MST
        mst_wt += blue_edge.first.first - red_edge.first.first;
        adjMatrix[red_edge.second.first][red_edge.second.second] = {-1, -1};
        adjMatrix[red_edge.second.second][red_edge.second.first] = {-1, -1};
        adjMatrix[blue_edge.second.first][blue_edge.second.second] = {blue_edge.first.first, 0};
        adjMatrix[blue_edge.second.second][blue_edge.second.first] = {blue_edge.first.first, 0};

        // Update parent and depth
        tie(parent, depth) = dfs(adjMatrix, 0);

        // Update red edge count
        red_edge_count -= (red_edge.first.second == 1);
    }

    cout << red_edge_count << endl << mst_wt << endl;
    return 0;
}