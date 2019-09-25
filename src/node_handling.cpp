#include "mapper.hpp"

/**
 * Creates a node based on parameters
 */
node create_node(const int cost, const int nswaps, const vector<edge> swaps) {
	node n;
	n.cost_fixed = cost;
	n.cost_heur  = n.lookahead_penalty = 0;
	n.total_cost = 0;
	n.qubits     = new int[positions];
	n.locations  = new int[nqubits];
	n.depths     = new int[positions];
	n.fidelities = new int[positions];
    n.nswaps     = nswaps;
	n.swaps      = swaps;
	n.done       = 1;
    return n;
}

/**
 * Creates a node
 */
node create_node() {
    return create_node(0, 0, vector<edge>());
}

/**
 * Creates a node based on a base node
 */
node create_node(const node& base, const edge& e) {
    node n = create_node(base.cost_fixed + COST_SWAP, base.nswaps, base.swaps);
    
	memcpy(n.qubits,     base.qubits,     sizeof(int) * positions);
	memcpy(n.locations,  base.locations,  sizeof(int) * nqubits);
	memcpy(n.depths,     base.depths,     sizeof(int) * positions);
    memcpy(n.fidelities, base.fidelities, sizeof(int) * positions);

	add_swap(n, e);

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

/*
 * adds a swap to the node
 */ 
void add_swap(node& n, const edge& e) {
	n.swaps.push_back(e);
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
	//clean up
	nodes.delete_queue();
	nodes = unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less>();
}
