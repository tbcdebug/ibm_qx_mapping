#include "mapper.hpp"

/**
 * Initializes the layer based on the gates of a circuit
 */
#if ONE_GATE_PER_LAYER
// this function generates one gate per layer
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
// this function generates the layers in a greedy fashion
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

/**
 * returns the index of the next layer with cnot gates
 */
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

/**
 * calculates the maximal width considering all layerss
 */
unsigned int calculate_max_layer_width() {
    unsigned int width = 0;
    for (vector<vector<QASMparser::gate>>::iterator it = layers.begin(); it != layers.end(); it++) {
		if (it->size() > width) {
			width = it->size();
		}
	}
    return width;
}
