#ifndef MAPPER_H_
#define MAPPER_H_


/**
 * Own Includes
 */
#include "qasm/QASMparser.h"
#include "unique_priority_queue.h"


/**
 * General Defines
 */
#define SUCCESS 0
#define ERROR   1


#define LOOK_AHEAD 1
#define HEURISTIC_ADMISSIBLE 0
#define USE_INITIAL_MAPPING 0 
#define MINIMAL_OUTPUT 1
#define DUMP_MAPPED_CIRCUIT 0

#define ARCH_LINEAR_N 0
#define ARCH_IBM_QX5 1

#ifndef ARCH
// assume default architecture
#define ARCH ARCH_LINEAR_N
#endif

/*
 * Constants
 */
// cost
const int COST_GATE     = 1;
const int COST_SWAP     = 7 * COST_GATE;





using namespace std;

extern double** dist;
extern int positions;
extern unsigned long ngates;
extern unsigned int nqubits;

struct edge {
	int v1;
	int v2;
};

inline bool operator<(const edge& lhs, const edge& rhs) {
	if (lhs.v1 != rhs.v1) {
		return lhs.v1 < rhs.v1;
	}
	return lhs.v2 < rhs.v2;
}

struct node {
	int cost_fixed;
	int cost_heur;
	int cost_heur2;
	int depth;
	int* qubits; // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int* locations; // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
	int nswaps;
	int done;
	vector<vector<edge>> swaps;
};

struct node_func_less {
	// true iff x < y
	bool operator()(const node& x, const node& y) const {
		for(int i=0; i < positions; i++) {
			if (x.qubits[i] != y.qubits[i]) {
				return x.qubits[i] < y.qubits[i];
			}
		}
		return false;
	}
};

struct node_cost_greater {
	// true iff x > y
	bool operator()(const node& x, const node& y) const {
		if ((x.cost_fixed + x.cost_heur + x.cost_heur2) != (y.cost_fixed + y.cost_heur + y.cost_heur2)) {
			return (x.cost_fixed + x.cost_heur + x.cost_heur2) > (y.cost_fixed + y.cost_heur + y.cost_heur2);
		}

		if(x.done == 1) {
			return false;
		}
		if(y.done == 1) {
			return true;
		}

		if (x.cost_heur + x.cost_heur2 != y.cost_heur + y.cost_heur2) {
			return x.cost_heur + x.cost_heur2 > y.cost_heur + y.cost_heur2;
		} else {
			return node_func_less{}(x, y);
		}

	}
};

struct cleanup_node {
	void operator()(const node& x) {
		delete[] x.qubits;
		delete[] x.locations;
	}
};

// circuit properties
struct circuit_properties {
	int* locations;
	int* qubits;
	int* depths;
	int* fidelities;
};

// dijkstra
struct dijkstra_node {
	int  pos;
	bool contains_correct_edge;
    int  length;
};

struct dijkstra_node_cmp {
	bool operator()(dijkstra_node* x, dijkstra_node* y) const {
		if(x->length != y->length) {
			return x->length > y->length;
		}

		if(!x->contains_correct_edge) {
			return true;
		}

		return y->contains_correct_edge;
	}
};

extern set<edge> graph;
extern vector<vector<QASMparser::gate> > layers;
extern unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less> nodes;

// coupling_graph
bool generate_graph(const string input);

// cost
double calculate_heuristic_cost(const dijkstra_node* node);

// layer_handling
vector<vector<QASMparser::gate>> init_layers(const vector<QASMparser::gate> &gates);
unsigned int get_next_layer(const unsigned int layer);
unsigned int calculate_max_layer_width();

// circuit_property_handling
circuit_properties create_circuit_properties();
void               delete_circuit_properties(circuit_properties& p);
void               adapt_circuit_properties(circuit_properties& p, const node& n);
void 			   update_properties(const int layer, circuit_properties& p);


int mapper(const vector<QASMparser::gate>& gates, vector<vector<QASMparser::gate>>& mapped_circuit, 
			vector<QASMparser::gate>& all_gates, int &total_swaps, circuit_properties& properties);

#endif /* MAPPER_H_ */