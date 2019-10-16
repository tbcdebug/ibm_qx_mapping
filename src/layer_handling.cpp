#include "mapper.hpp"

#if ONE_GATE_PER_LAYER
vector<vector<QASMparser::gate>> init_layers(const vector<QASMparser::gate> &gates) {
	unsigned int layer = 0;
	vector<vector<QASMparser::gate>> layer_gates;
	for (vector<QASMparser::gate>::const_iterator it = gates.begin(); it != gates.end(); it++) {
		layer_gates.push_back(vector<QASMparser::gate>());
		layer_gates[layer].push_back(*it);	
		layer++;		
	}
	return layer_gates;
}
#else
vector<vector<QASMparser::gate>> init_layers(const vector<QASMparser::gate> &gates) {
	vector<vector<QASMparser::gate>> layer_gates;
	unsigned int layer;
	int* last_layer = new int[positions];
	for (int i = 0; i < positions; i++) {
		last_layer[i] = -1;
	}

	for (vector<QASMparser::gate>::const_iterator it = gates.begin(); it != gates.end(); it++) {
		QASMparser::gate g = *it;
		if(g.control == -1) { // single qubit operation
			layer                = last_layer[g.target] + 1;
			last_layer[g.target] = layer;                                	
		} else {			  // control also available
			layer                = max(last_layer[g.target], last_layer[g.control]) + 1; 
			last_layer[g.target] = last_layer[g.control] = layer;// set layer for both operands
		}

		if (layer_gates.size() <= layer) {									// append new layer if not enough layers
			layer_gates.push_back(vector<QASMparser::gate>());
		}
		layer_gates[layer].push_back(g);			
	}
	delete[] last_layer;
	return layer_gates;
}
#endif

// returns the next layer with cnot gates
unsigned int get_next_layer(const unsigned int layer) {
	unsigned int next_layer = layer+1;
	while(next_layer < layers.size()) {
		for(vector<QASMparser::gate>::iterator it = layers[next_layer].begin(); it != layers[next_layer].end(); it++) {
			if(it->control != -1) {
				return next_layer;
			}
		}
		next_layer++;
	}
	return -1;
}

/*
void map_unmapped_gates(const int layer, circuit_properties& p, node& n, vector<int>& considered_qubits) {
	int* loc = p.locations;
	int* map = p.qubits;
	
	for (vector<QASMparser::gate>::iterator it = layers[layer].begin(); it != layers[layer].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			considered_qubits.push_back(g.control);
			considered_qubits.push_back(g.target);
			if(loc[g.control] == -1 && loc[g.target] == -1) { //TODO better mapping
				set<edge> possible_edges;
				for(set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
					if(map[it->v1] == -1 && map[it->v2] == -1) {
						possible_edges.insert(*it);
					}
				}
				if(!possible_edges.empty()) {
					edge e         = *possible_edges.begin();
					loc[g.control] = e.v1;
					map[e.v1]      = g.control;
					loc[g.target]  = e.v2;
					map[e.v2]      = g.target;
				} else {
					cerr << "no edge available!";
					exit(1);
				}
			} else if(loc[g.control] == -1) { //DONE min function
				map_to_min_distance(map, loc, g.target, g.control);
			} else if(loc[g.target] == -1) {
				map_to_min_distance(map, loc, g.control, g.target);
			}
			n.cost_heur = max(n.cost_heur, dist[loc[g.control]][loc[g.target]]);
		}
	}
}
*/
// calculates the width of the layer
unsigned int calculate_max_layer_width() {
    unsigned int width = 0;
    for (vector<vector<QASMparser::gate>>::iterator it = layers.begin(); it != layers.end(); it++) {
		if (it->size() > width) {
			width = it->size();
		}
	}
    return width;
}
