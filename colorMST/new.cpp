
#include<bits/stdc++.h>
  
using namespace std; 
vector<vector<pair<int, int>>> adjMatrix;
// Creating shortcut for an integer pair 
typedef pair<int, int> iPair;
 
// Structure to represent a graph 
struct Graph 
{ 
	int V, E;
 	vector< pair<pair<int,int>, iPair> > edges; 
	// Constructor 
	Graph(int V, int E) 
	{ 
		this->V = V; 
		this->E = E; 
 	} 

	// Utility function to add an edge  // red=1 means red
	void addEdge(int u, int v, int w, int red) 
	{ 
		edges.push_back({{w,red}, {u, v}});
	} 

}; 
// To represent Disjoint Sets 
struct DisjointSets 
{ 
	int *parent, *rnk; 
	int n; 

	// Constructor. 
	DisjointSets(int n) 
	{ 
		// Allocate memory 
		this->n = n; 
		parent = new int[n+1]; 
		rnk = new int[n+1]; 

		// Initially, all vertices are in 
		// different sets and have rank 0. 
		for (int i = 0; i <= n; i++) 
		{ 
			rnk[i] = 0; 

			//every element is parent of itself 
			parent[i] = i; 
		} 
	} 

	// Find the parent of a node 'u' 
	// Path Compression 
	int find(int u) 
	{ 
		/* Make the parent of the nodes in the path 
		from u--> parent[u] point to parent[u] */
		if (u != parent[u]) 
			parent[u] = find(parent[u]); 
		return parent[u]; 
	} 

	// Union by rank 
	void merge(int x, int y) 
	{ 
		x = find(x), y = find(y); 

		/* Make tree with smaller height 
		a subtree of the other tree */
		if (rnk[x] > rnk[y]) 
			parent[y] = x; 
		else // If rnk[x] <= rnk[y] 
			parent[x] = y; 

		if (rnk[x] == rnk[y]) 
			rnk[y]++; 
	} 
}; 
// Function to perform Kruskal's algorithm and return MST edges and weight

// Function to perform Kruskal's algorithm and return MST edges and weight
pair<int, vector<pair<pair<int, int>, iPair>>> kruskalMST(Graph& g) 
{ 
	int mst_wt = 0; // Initialize result 
	vector<pair<pair<int, int>, iPair>> mst_edges; // To store edges of MST

	// Sort edges in increasing order on basis of cost 
	sort(g.edges.begin(), g.edges.end()); 

	// Create disjoint sets 
	DisjointSets ds(g.V); 

	// Iterate through all sorted edges 
	for (auto it = g.edges.begin(); it != g.edges.end(); it++) 
	{ 
		int u = it->second.first; 
		int v = it->second.second; 

		int set_u = ds.find(u); 
		int set_v = ds.find(v); 

		// Check if the selected edge is creating 
		// a cycle or not (Cycle is created if u 
		// and v belong to the same set) 
		if (set_u != set_v) 
		{ 
			// Current edge will be in the MST 
			mst_edges.push_back(*it); 

			// Update MST weight 
			mst_wt += it->first.first; 

			// Merge two sets 
			ds.merge(set_u, set_v); 
		} 
	} 

	return {mst_wt, mst_edges}; 
}


pair<vector<int>, vector<int>> dfs(int startVertex) {
	int V = adjMatrix.size();
	vector<int> parent(V, -1);
	vector<int> depth(V, -1);
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

vector<pair<pair<int,int>,pair<int,int>>> retruncycleedges(pair<vector<int>, vector<int>>parentdepth,const pair<pair<int, int>, iPair>& new_edge){
int vertex1=new_edge.second.first;
int vertex2=new_edge.second.second;
auto parent=parentdepth.first;
auto depth=parentdepth.second;
vector<pair<pair<int,int>,pair<int,int>>>cycleedges;
while(depth[vertex1]>depth[vertex2]){
cycleedges.push_back({adjMatrix[vertex1][parent[vertex1]],{vertex1,parent[vertex1]}});
	vertex1=parent[vertex1];

}
while(depth[vertex2]>depth[vertex1]){
	cycleedges.push_back({adjMatrix[vertex2][parent[vertex2]],{vertex2,parent[vertex2]}});
	vertex2=parent[vertex2];
	}
while(vertex1!=vertex2){
cycleedges.push_back({adjMatrix[vertex1][parent[vertex1]],{vertex1,parent[vertex1]}});
cycleedges.push_back({adjMatrix[vertex2][parent[vertex2]],{vertex2,parent[vertex2]}});
vertex1=parent[vertex1];
vertex2=parent[vertex2];
}
return cycleedges;
	}


pair<iPair, iPair> findMaxWeightRedEdgeInCycle(const vector<pair<pair<int, int>, pair<int, int>>>& cycle_edges) {
		pair<iPair, iPair> max_red_edge = {{-1, 1}, {-1, -1}}; // Default value if no red edge is found
		for (const auto& edge : cycle_edges) {
			if (edge.first.second == 1) { // Check if the edge is red
				if (edge.first.first > max_red_edge.first.first) {
					max_red_edge = edge;
				}
			}
		}
		return max_red_edge;
	}

	pair<iPair, iPair> processNewEdgeWithAdjMatrix( pair<vector<int>, vector<int>> parentdepth, const pair<pair<int, int>, iPair>& new_edge) {
		// Find the cycle edges with the new edge added
		// Check if the edge is already in the adjacency matrix
		int u = new_edge.second.first;
		int v = new_edge.second.second;
		if (adjMatrix[u][v].first != -1) {
			return {{-1, -1}, {-1, -1}}; // Return default value if edge is already present
		}
		auto cycle_edges = retruncycleedges( parentdepth, new_edge);

		// Find the maximum weight red edge in the cycle
		return findMaxWeightRedEdgeInCycle(cycle_edges);
	}




int main() 
{ 
	int V, E; 
	int threshold;  

	cin >> V;
	cin >> E;
	cin >> threshold;
	Graph g(V, E); 
 

	int u, v, w, r;

	for (int i=0; i< E; i++){
		cin >> u;
		cin >> v;
		cin >> w;
		cin >> r;
		g.addEdge(u, v, w, r); 
	}
	auto result = kruskalMST(g);

	int red_edge_count = 0;
	for (const auto& edge : result.second) {
		if (edge.first.second == 1) {
			red_edge_count++;
		}
	}
	// Create an adjacency matrix for the MST
	
	adjMatrix.resize(V, vector<pair<int, int>>(V, {-1, -1}));
	for (const auto& edge : result.second) {
		int u = edge.second.first;
		int v = edge.second.second;
		int weight = edge.first.first;
		int isRed = edge.first.second;
		adjMatrix[u][v] = {weight, isRed};
		adjMatrix[v][u] = {weight, isRed}; // Since the graph is undirected
	}

	// Print the adjacency matrix

	vector<pair<pair<int, int>, iPair>> blue_edges_not_in_mst;
	for (const auto& edge : g.edges) {
		if (edge.first.second == 0) { // Check if the edge is blue
			bool in_mst = false;
			for (const auto& mst_edge : result.second) {
				if (edge == mst_edge) {
					in_mst = true;
					break;
				}
			}
			if (!in_mst) {
				blue_edges_not_in_mst.push_back(edge);
			}
		}
	}
	// Sort blue edges not in MST by weight
	sort(blue_edges_not_in_mst.begin(), blue_edges_not_in_mst.end());
	
	
	




bool xoxoxoxo=true;
auto parentdepth=dfs(0);
	for(int pp=0;pp<=blue_edges_not_in_mst.size();pp++)	{vector<pair<pair<iPair, iPair>,pair<iPair, iPair>>>blueredmap;
		
		if(xoxoxoxo) {parentdepth=dfs(0);xoxoxoxo=false;}

		for (const auto& blue_edge : blue_edges_not_in_mst) {
		
			auto max_red_edge = processNewEdgeWithAdjMatrix( parentdepth, blue_edge);
			
					
			
			if (max_red_edge.second.first != -1 && max_red_edge.second.second != -1) {
				
				int new_mst_weight = result.first + blue_edge.first.first - max_red_edge.first.first;
				
				if (new_mst_weight <= threshold) {	
					
				blueredmap.push_back({blue_edge, max_red_edge});
			}}
			
		}
	
	if(blueredmap.size()!=0){
	int modifind=0;
	int xoxo=blueredmap[0].first.first.first-blueredmap[0].second.first.first;
	for (int i=0; i<blueredmap.size(); i++){
		auto bluered=blueredmap[i];
		if(xoxo>bluered.first.first.first-bluered.second.first.first){
			xoxo=bluered.first.first.first-bluered.second.first.first;
			modifind=i;
		}
	}
	auto bluered=blueredmap[modifind];
	// Replace the red edge with the blue edge in the MST
	bool xo=false;
	

			adjMatrix[bluered.second.second.first][bluered.second.second.second] = {-1, -1};
			adjMatrix[bluered.second.second.second][bluered.second.second.first] = {-1, -1};
			
		
	
	
	
	result.first=result.first+xoxo;
	adjMatrix[bluered.first.second.first][bluered.first.second.second] = {bluered.first.first.first, 0};
	adjMatrix[bluered.first.second.second][bluered.first.second.first] = {bluered.first.first.first, 0};
xoxoxoxo=true;	
}
else{
	break;}



	}





	// Recalculate the MST weight and red edges from the adjacency matrix
	result.first = 0;
	result.second.clear();
	red_edge_count = 0;

	for (int i = 0; i < V; ++i) {
		for (int j = i + 1; j < V; ++j) {
			if (adjMatrix[i][j].first != -1) {
				result.first += adjMatrix[i][j].first;
				result.second.push_back({{adjMatrix[i][j].first, adjMatrix[i][j].second}, {i, j}});
				if (adjMatrix[i][j].second == 1) {
					red_edge_count++;
				}
			}
		}
	}

	// Print the final MST weight and number of red edges
	
	
	cout << red_edge_count << endl;
	cout <<result.first << endl;
 

	return 0; 
}


