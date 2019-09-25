#include "mapper.hpp"

#include "math.h"


int get_maximal_depth(const int* depths) {	
	// calculate max 
    int max_depth = 0;
	for(int i = 0; i < positions; i++) {
		max_depth = max(max_depth, depths[i]);
	}
    return max_depth;
}

long long fidelity_cost(const int* fidelities) {
	int       avg              = 0;
	int       nqubits_not_null = 0;
	long long variance         = 0;
	for(int i = 0; i < positions; i++) {  // calcualte average
		if(fidelities[i] != 0) {
			avg              += fidelities[i];
			nqubits_not_null ++;
		}
	}
	avg /= nqubits_not_null;
	for(int i = 0; i < positions; i++) {
		if(fidelities[i] != 0) {
			long long diff = fidelities[i] - avg;
			//variance += abs(diff); 
			variance += diff * diff;
		}
	}
	
	return sqrt(variance/nqubits_not_null);
}

double get_total_cost(const node& n) {
#if SPECIAL_OPT
	return (fidelity_cost(n.fidelities)                       * FIDELITY_NORM)    + 
		   (get_maximal_depth(n.depths)/((double) DEPTH_SWAP) * DEPTH_PERCENTAGE) + 
		   (n.cost_fixed/((double)  COST_SWAP)                * COST_PERCENTAGE);
#else
    return n.cost_fixed;
#endif
}

double calculate_heuristic_cost(const dijkstra_node* node) {
	int path_length = node->length - 1;
	if(node->contains_correct_edge) {
#if SPECIAL_OPT
		return path_length;
#else 
		return path_length * COST_SWAP;
#endif
	}
#if SPECIAL_OPT
	return path_length + INVERSE; 
#else
	return path_length * COST_SWAP + 4;
#endif
}


double heuristic_function(const double old_heur, const double new_heur) {
#if HEURISTIC_ADMISSIBLE
	return max(old_heur, new_heur);
#else
	return old_heur + new_heur;
#endif
}

double get_heuristic_cost(const double cost_heur, const node& n, const QASMparser::gate& g) {
	return heuristic_function(cost_heur, dist[n.locations[g.control]][n.locations[g.target]]);
}