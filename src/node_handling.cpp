#include "mapper.hpp"

#include <cstring>

static void apply_edge(const node& n, const edge e) {
	int tmp_qubit1 = n.qubits[e.v1];
	int tmp_qubit2 = n.qubits[e.v2];

	n.qubits[e.v1] = tmp_qubit2;
	n.qubits[e.v2] = tmp_qubit1;

	if (tmp_qubit1 != -1) {
		n.locations[tmp_qubit1] = e.v2;
	}
	if (tmp_qubit2 != -1) {
		n.locations[tmp_qubit2] = e.v1;
	}

	int max_depth = max(n.depths[e.v1], n.depths[e.v2]) + DEPTH_SWAP;
	n.depths[e.v1]      = max_depth;
	n.fidelities[e.v1] += FIDELITY_SWAP;
	n.depths[e.v2]      = max_depth;
	n.fidelities[e.v2] += FIDELITY_SWAP;
}

/**
 * Creates a node based on parameters
 */
node create_node(const int cost, const int nswaps, SWAP_LIST_TYPE swaps) {
	node n;
	n.cost_fixed = cost;
	n.cost_heur  = n.lookahead_penalty = 0;
	n.total_cost = 0;
	n.qubits     = new int[positions];
	n.locations  = new int[nqubits];
	n.depths     = new int[positions];
	n.fidelities = new int[positions];
    n.nswaps     = nswaps;
	n.done       = 1;
    n.swaps      = SWAP_LIST_TYPE();
	for (SWAP_LIST_TYPE::iterator it = swaps.begin(); it != swaps.end(); it++) {
		n.swaps.push_back(*it);
	}
	
    return n;
}

/**
 * Creates a node
 */
node create_node() {
    return create_node(0, 0, SWAP_LIST_TYPE());
}

/**
 * Creates a node based on a base node
 */
node create_node(const node& base, const edge* new_swaps, const int nswaps) {
    node n = create_node(base.cost_fixed + COST_SWAP * nswaps, base.nswaps + nswaps, base.swaps);
    
	memcpy(n.qubits,     base.qubits,     sizeof(int) * positions);
	memcpy(n.locations,  base.locations,  sizeof(int) * nqubits);
	memcpy(n.depths,     base.depths,     sizeof(int) * positions);
    memcpy(n.fidelities, base.fidelities, sizeof(int) * positions);


    SWAP_TYPE n_swaps;

	for (int i = 0; i < nswaps; i++) {
		apply_edge(n, new_swaps[i]);
        n_swaps.push_back(new_swaps[i]);
	}
    
    n.swaps.push_back(n_swaps);

    /*
    	
	n.swaps.push_back(e);
    */

	n.total_cost = get_total_cost(n);

    return n;
}

/*
 * updates the node based on the circuit properties
 */
void update_node(node& n, const circuit_properties& p) {
	memcpy(n.qubits,     p.qubits,     sizeof(int) * positions);
	memcpy(n.locations,  p.locations,  sizeof(int) * nqubits);
	memcpy(n.depths,     p.depths,     sizeof(int) * positions);
	memcpy(n.fidelities, p.fidelities, sizeof(int) * positions);
}

/**
 * Checks if a node is a goal and stops
 */
void check_if_not_done(node& n, const int value) {
#if SPECIAL_OPT
	if(value >= 1) {
#else
	if(value > 4) {
#endif
		n.done = 0;
	}
}

/**
 * Deletes a node
 */
void delete_node(const node& n) {
   cleanup_node()(n);
}

/**
 * Deletes all nodes
 */
void delete_nodes() {
    nodes.delete_queue();
    nodes = unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less>();
}
