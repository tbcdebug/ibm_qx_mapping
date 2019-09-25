#ifndef MAPPER_H_
#define MAPPER_H_

#include <iostream>
#include <algorithm>
#include <string.h>
#include <set>
#include <climits>
#include <fstream>

#include "qasm/QASMparser.h"
#include "unique_priority_queue.h"


/**
 * General Defines
 */
#define SUCCESS 0
#define ERROR   1


/**
 * Control Defines
 */
#define LOOK_AHEAD 1
#define HEURISTIC_ADMISSIBLE 0
#define USE_INITIAL_MAPPING 1

#define ARCH_LINEAR_N 0
#define ARCH_IBM_QX5 1

/*
#ifndef ARCH
// assume default architecture
#define ARCH ARCH_LINEAR_N
#endif
*/

using namespace std;

/*
 * Constants
 */
// fidelity
const int FIDELITY_GATE = 1;
const int FIDELITY_CNOT = 5;
const int FIDELITY_SWAP = 2 * FIDELITY_GATE + 3 * FIDELITY_CNOT;

// depth
const int DEPTH_GATE    = 1;
const int DEPTH_SWAP    = 5 * DEPTH_GATE;

// cost
const int COST_GATE     = 1;
const int COST_SWAP     = 7 * COST_GATE;

// cost factors
const double COST_PERCENTAGE  = 1;//0.05;
const double DEPTH_PERCENTAGE = 1 - COST_PERCENTAGE;
const double FIDELITY_FACTOR  = 0;
const double FIDELITY_NORM    = FIDELITY_FACTOR / 1000;

const double INVERSE = DEPTH_PERCENTAGE * (((double)2) * DEPTH_GATE / DEPTH_SWAP) + COST_PERCENTAGE  * 0.57;

// lookahead
const int    N_LOOK_AHEADS             = 1;
const double FIRST_LOOK_AHEAD_FACTOR   = 0.9;
const double GENERAL_LOOK_AHEAD_FACTOR = 0.5;

// nn positions - only used when nn coupling graph should be used
const int NN_POSITIONS = 16;


/*
 * GLOBAL Variables
 */
extern double**      dist;
extern int           positions;
extern unsigned long ngates;
extern unsigned int  nqubits;
 
extern unsigned int  max_node_size;

/*
 * Graph Struct
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

/* 
 * Node Handling
 */
struct node {
	int    cost_fixed;
	double cost_heur;
	double lookahead_penalty; // look ahead penalty
	double total_cost;         
	int*   qubits;             // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int*   locations;          // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
	int*   depths;
	int*   fidelities;
	int    nswaps;
	int    done;
	vector<edge> swaps;
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

extern set<edge>                                                                    graph;
extern vector<vector<QASMparser::gate>>                                             layers;
extern unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less> nodes;

/*
 * Function Declarations
 */

// utility
bool contains(const vector<int>& v, const int e);
void map_to_min_distance(int* map, int* loc, const int source, const int target);
void generate_circuit(vector<vector<QASMparser::gate>>& mapped_circuit, const vector<QASMparser::gate>& all_gates);
void fix_positions_of_single_qubit_gates(int* locations, int* qubits, vector<QASMparser::gate>& all_gates);

// coupling_graph
void   set_dijkstra_node(dijkstra_node* nodes, priority_queue<dijkstra_node*, vector<dijkstra_node*>, dijkstra_node_cmp>& queue,
					   const int parent, const int pos, const bool contains_correct_edge);
void   dijkstra(dijkstra_node* nodes, const int start, const set<edge>& graph);
void   build_dist_table(const set<edge>& graph);
bool   generate_graph(const string input);

void   build_graph_linear(int nqubits);
void   build_graph_NN(int nqubits);
void   build_graph_QX5();

// node_handling
node create_node();
node create_node(const node& base, const edge& e);
void update_node(node& n, const circuit_properties& p);
void add_swap(node& n, const edge& e);
void check_if_not_done(node& n, const int value);
void delete_node(const node& n);
void delete_nodes(); 

// mapping
void lookahead(const int layer, node& new_node);
void expand_node_add_one_swap(const vector<int>& qubits, const edge e, node base_node, 
							  const int layer);
void expand_node(const vector<int>& qubits,	const node base_node, const int layer);
node a_star_fixlayer(const int layer, circuit_properties &properties);
void mapper(const vector<QASMparser::gate>& gates, vector<vector<QASMparser::gate>>& mapped_circuit, 
			vector<QASMparser::gate>& all_gates, int &total_swaps, circuit_properties& properties);

// layer_handling
vector<vector<QASMparser::gate>> init_layers(const vector<QASMparser::gate>& gates);
unsigned int                     get_next_layer(const unsigned int layer);
void                             map_unmapped_gates(const int layer, circuit_properties& p, 
                                                    node& n, vector<int>& considered_qubits);
unsigned int                     calculate_layer_width();

// circuit_property_handling
circuit_properties create_circuit_properties();
void               delete_circuit_properties(circuit_properties& p);
void               adapt_circuit_properties(circuit_properties& p, const node& n);
void 			   update_properties(const int layer, circuit_properties& p);

// cost
int       get_maximal_depth(const int* depths);
int       depth_cost(const node& n);
long long fidelity_cost(const int* fidelities);
double    get_total_cost(const node& n);
double    calculate_heuristic_cost(const dijkstra_node* node);
double    heuristic_function(const double old_heur, const double new_heur);
double    get_heuristic_cost(const double cost_heur, const node& n, const QASMparser::gate& g);

#endif /* MAPPER_H_ */