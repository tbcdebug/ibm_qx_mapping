#ifndef MAPPER_H_
#define MAPPER_H_


#include <vector>


using namespace std;


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

#define ARCH_LINEAR_N 0
#define ARCH_IBM_QX5 1
#define ARCH_LINEAR_NN 2

#define UNUSED_FUNCTION __attribute__ ((unused))

/**
 * Control Defines
 */

#define LOOK_AHEAD           1
#define HEURISTIC_ADMISSIBLE 0
#define USE_INITIAL_MAPPING  0
#define MINIMAL_OUTPUT       1
#define DUMP_MAPPED_CIRCUIT  0

#ifndef ARCH
#define ARCH ARCH_LINEAR_N     // assume default architectures
#endif

/*
 * Constants
 */
// cost
const int COST_GATE     = 1;
const int COST_SWAP     = 7 * COST_GATE;

// fidelity
const int FIDELITY_GATE = 1;
const int FIDELITY_CNOT = 5;
const int FIDELITY_SWAP = 2 * FIDELITY_GATE + 3 * FIDELITY_CNOT;

// depth
const int DEPTH_GATE    = 1;
const int DEPTH_SWAP    = 5 * DEPTH_GATE;





extern double**      dist;
extern int           positions;
extern unsigned long ngates;
extern unsigned int  nqubits;

/**
 * Types and structs
 */

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

typedef vector<edge>      SWAP_TYPE;
typedef vector<SWAP_TYPE> SWAP_LIST_TYPE;   

struct node {
	int    cost_fixed;
	double cost_heur;
	double lookahead_penalty;
	double total_cost;
	int    depth;
	int*   qubits;    // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int*   locations; // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
	int*   depths;
	int*   fidelities;
	int    nswaps;
	int    done;
	SWAP_LIST_TYPE swaps;
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
		if ((x.cost_fixed + x.cost_heur + x.lookahead_penalty) != (y.cost_fixed + y.cost_heur + y.lookahead_penalty)) {
			return (x.cost_fixed + x.cost_heur + x.lookahead_penalty) > (y.cost_fixed + y.cost_heur + y.lookahead_penalty);
		}

		if(x.done == 1) {
			return false;
		}
		if(y.done == 1) {
			return true;
		}

		if (x.cost_heur + x.lookahead_penalty != y.cost_heur + y.lookahead_penalty) {
			return x.cost_heur + x.lookahead_penalty > y.cost_heur + y.lookahead_penalty;
		} else {
			return node_func_less{}(x, y);
		}

	}
};

struct cleanup_node {
	void operator()(const node& n) {
		delete[] n.locations;
		delete[] n.qubits;
		delete[] n.depths;
		delete[] n.fidelities;
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
double get_total_cost(const node& n);
double heuristic_function(const double old_heur, const double new_heur);
double get_heuristic_cost(const double cost_heur, const node& n, const QASMparser::gate& g);

// node_handling
node create_node();
node create_node(const node& base, const edge* new_swaps, const int nswaps);
void update_node(node& n, const circuit_properties& p);
void add_swap(node& n, const edge& e);
void check_if_not_done(node& n, const int value);
void delete_node(const node& n);
void delete_nodes(); 

// layer_handling
vector<vector<QASMparser::gate>> init_layers(const vector<QASMparser::gate> &gates);
unsigned int get_next_layer(const unsigned int layer);
unsigned int calculate_max_layer_width();

// circuit_property_handling
circuit_properties create_circuit_properties();
void               delete_circuit_properties(circuit_properties& p);
void               adapt_circuit_properties(circuit_properties& p, const node& n);
void 			   update_properties(const int layer, circuit_properties& p);


void mapping(const vector<QASMparser::gate>& gates, vector<vector<QASMparser::gate>>& mapped_circuit, 
			vector<QASMparser::gate>& all_gates, int &total_swaps, circuit_properties& properties);

#endif /* MAPPER_H_ */