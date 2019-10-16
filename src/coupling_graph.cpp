#include "mapper.hpp"

static void set_dijkstra_node(dijkstra_node* nodes, priority_queue<dijkstra_node*, vector<dijkstra_node*>, dijkstra_node_cmp>& queue,
			       		      const int parent, const int pos, const bool contains_correct_edge);
static void dijkstra(dijkstra_node* nodes, const int start, const set<edge>& graph);
static void build_dist_table(const set<edge>& graph);
static void build_graph_linear(int nqubits)           UNUSED_FUNCTION;
static void build_graph_NN(int nqubits)               UNUSED_FUNCTION;
static void build_graph_QX5()                         UNUSED_FUNCTION;


static void set_dijkstra_node(dijkstra_node* nodes, priority_queue<dijkstra_node*, vector<dijkstra_node*>, dijkstra_node_cmp>& queue,
					          const int parent, const int pos, const bool contains_correct_edge) {
	if(nodes[pos].length < 0) {
		nodes[pos].contains_correct_edge = contains_correct_edge;
		nodes[pos].length                = nodes[parent].length + 1;
		queue.push(nodes + pos);
	} else if(!nodes[pos].contains_correct_edge && contains_correct_edge && nodes[pos].length == nodes[parent].length + 1) {
		priority_queue<dijkstra_node*, vector<dijkstra_node*>, dijkstra_node_cmp> temp;
		nodes[pos].contains_correct_edge = true;

		while(!queue.empty()) {
			temp.push(queue.top());
			queue.pop();
		}
		queue = temp;
	}
}

static void dijkstra(dijkstra_node* nodes, const int start, const set<edge>& graph) {
	priority_queue<dijkstra_node*, vector<dijkstra_node*>, dijkstra_node_cmp> queue;
	queue.push(nodes + start);

	while(!queue.empty()) {
		dijkstra_node* current = queue.top();
		queue.pop();
		int cur = current->pos;
		for (set<edge>::iterator it = graph.begin(); it != graph.end();
				it++) {
			edge e = *it;
			if (cur == e.v1) { 
				set_dijkstra_node(nodes, queue, e.v1, e.v2, true);	
			} else if (cur == e.v2) {
				set_dijkstra_node(nodes, queue, e.v2, e.v1, current->contains_correct_edge);
			}
		}
	}
}

// build the distance table
static void build_dist_table(const set<edge>& graph) {
	dist = new double*[positions];

	for (int i = 0; i < positions; i++) {
		dist[i] = new double[positions];
	}

	for (int i = 0; i < positions; i++) {
		dijkstra_node* nodes = new dijkstra_node[positions];
		for (int j = 0; j < positions; j++) {
			nodes[j].pos                   = j;
			nodes[j].contains_correct_edge = false;
			nodes[j].length                = -1;
		}
		nodes[i].length = 0;
		
		dijkstra(nodes, i, graph);
		for (int j = 0; j < positions; j++) {
			if (i != j) {
				dist[i][j] = calculate_heuristic_cost(nodes + j);
			} else {
				dist[i][j] = 0;
			}
		}

		delete[] nodes;
	}
}

// generate graph from the file, if the file is empty a QX5 graph is generated
bool generate_graph(const string input) {
    if(input.empty()) {
#if ARCH == ARCH_LINEAR_N
		build_graph_linear(nqubits);
#elif ARCH == ARCH_LINEAR_NN
		build_graph_NN(nqubits):
#elif ARCH == ARCH_IBM_QX5
		build_graph_QX5();
#else
    	static_assert(false, "No architecture specified!");
#endif
    } else {
        graph.clear();
		ifstream infile(input, ifstream::in);
		if(infile.fail()) {
			cerr << "Failed to open file '" << input << "'!" << endl;
			return false;
		} else {
			string line;
			if(getline(infile, line)) {
				sscanf(line.c_str(), "Positions: %d", &positions);
				cout << "Positions: " << positions << endl;
			} else {
				cout << "First Line has to be: Positions: [0-9]*" << endl;
				return false;
			}
			while (getline(infile, line)) {
				int e1, e2;
				edge e;
				if (sscanf(line.c_str(), "[%d,%d]", &e1, &e2) == 2) {
					//cout << "Edge: " << e1 << "  " << e2 << "  detected." << endl;
					e.v1 = e1;
					e.v2 = e2;
					graph.insert(e);
				}
			}
			infile.close();
		}	
		cout << "Finished reading of the coupling graph" << endl;
    }

    build_dist_table(graph);
	cout << "Finished building of the dist table" << endl;
	return true;
}

//build a graph representing the coupling map of nearest neighbor
static void build_graph_linear(int nqubits) {
	graph.clear();
    
	positions = nqubits;
	for(int i = 0; i < nqubits - 1; i++) {
        graph.insert(edge{i, i+1});
        graph.insert(edge{i+1, i});
	}
}

//build a graph representing the coupling map of IBM QX5
static void build_graph_QX5() {
	graph.clear();
	positions = 16;
	edge e;
	e.v1 = 1;
	e.v2 = 0;
	graph.insert(e);
	e.v1 = 1;
	e.v2 = 2;
	graph.insert(e);
	e.v1 = 2;
	e.v2 = 3;
	graph.insert(e);
	e.v1 = 3;
	e.v2 = 14;
	graph.insert(e);
	e.v1 = 3;
	e.v2 = 4;
	graph.insert(e);
	e.v1 = 5;
	e.v2 = 4;
	graph.insert(e);
	e.v1 = 6;
	e.v2 = 5;
	graph.insert(e);
	e.v1 = 6;
	e.v2 = 11;
	graph.insert(e);
	e.v1 = 6;
	e.v2 = 7;
	graph.insert(e);
	e.v1 = 7;
	e.v2 = 10;
	graph.insert(e);
	e.v1 = 8;
	e.v2 = 7;
	graph.insert(e);
	e.v1 = 9;
	e.v2 = 8;
	graph.insert(e);
	e.v1 = 9;
	e.v2 = 10;
	graph.insert(e);
	e.v1 = 11;
	e.v2 = 10;
	graph.insert(e);
	e.v1 = 12;
	e.v2 = 5;
	graph.insert(e);
	e.v1 = 12;
	e.v2 = 11;
	graph.insert(e);
	e.v1 = 12;
	e.v2 = 13;
	graph.insert(e);
	e.v1 = 13;
	e.v2 = 4;
	graph.insert(e);
	e.v1 = 13;
	e.v2 = 14;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 0;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 14;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 2;
	graph.insert(e);
}

static void build_graph_NN(int nqubits) {
	build_graph_linear(nqubits);
	positions = 16;
}
