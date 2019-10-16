#include "mapper.hpp"

#include <string.h>

void initial_mapping(circuit_properties& properties) {
	int* qubits    = properties.qubits;
	int* locations = properties.locations;

	for (vector<QASMparser::gate>::iterator it = layers[0].begin(); it != layers[0].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			for(set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
				if(qubits[it->v1] == -1 && qubits[it->v2] == -1) {
					qubits[it->v1]       = g.control;
					qubits[it->v2]       = g.target;
					locations[g.control] = it->v1;
					locations[g.target]  = it->v2;
					break;
				}
			}
		}
	}
	for(unsigned int i = 0; i < nqubits; i++) {
		if(locations[i] == -1) {
			int j = 0;
			while(qubits[j] != -1){
				j++;
			}
			locations[i] = j;
			qubits[j]    = i;
		}
	}
}

// maps the qubit to the physical location that is not yet mapped and has minimum distance 
void map_to_min_distance(int* map, int* loc, const int source, const int target) {
	int min     = 1000;
	int min_pos = -1;
	for(int i = 0; i < positions; i++) {
		if(map[i] == -1 && dist[loc[source]][i] < min) {
			min     = dist[loc[source]][i];
			min_pos = i;
		}
	}
	map[min_pos] = target;
	loc[target]  = min_pos;
}

void map_unmapped_gates(const int layer, circuit_properties& p, node& n, vector<int>& considered_qubits) {
	int* map = p.qubits;
	int* loc = p.locations;

	//Find a mapping for all logical qubits in the CNOTs of the layer that are not yet mapped
	for (vector<QASMparser::gate>::iterator it = layers[layer].begin(); it != layers[layer].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			considered_qubits.push_back(g.control);
			considered_qubits.push_back(g.target);

			if(loc[g.control] == -1 && loc[g.target] == -1) {
                set<edge> possible_edges;
				for(set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
					if(map[it->v1] == -1 && map[it->v2] == -1) {
						possible_edges.insert(*it);
					}
				}
				if(!possible_edges.empty()) {
					edge e = *possible_edges.begin();
					loc[g.control] = e.v1;
					map[e.v1] = g.control;
					loc[g.target] = e.v2;
					map[e.v2] = g.target;
				} else {
                    cerr << "no edge available!";
                    exit(1);
				}
			} else if(loc[g.control] == -1) {
				map_to_min_distance(map, loc, g.target, g.control);
			} else if(loc[g.target] == -1) {
				map_to_min_distance(map, loc, g.control, g.target);
			}
			n.cost_heur = max(n.cost_heur, dist[loc[g.control]][loc[g.target]]);
		}
	}
}

// fix the position of the single qubit gates
void fix_positions_of_single_qubit_gates(int* locations, int* qubits, vector<QASMparser::gate>& all_gates) {
	for(vector<QASMparser::gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];

			qubits[it->control] = tmp_qubit2;
			qubits[it->target]  = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		}
		if(it->target < 0) {
			int target = -(it->target +1);
			it->target = locations[target];
			if(locations[target] == -1) {
				//This qubit occurs only in single qubit gates -> it can be mapped to an arbirary physical qubit
				int loc = 0;
				while(qubits[loc] != -1) {
					loc++;
				}
				locations[target] = loc;
			}
		}
	}
}

void generate_circuit(vector<vector<QASMparser::gate>>& mapped_circuit, const vector<QASMparser::gate>& all_gates) {
	int *last_layer = new int[positions];
	for(int i = 0; i < positions; i++) {
		last_layer[i] = -1;
	}

	//build resulting circuit
	for(vector<QASMparser::gate>::const_iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			continue;
		}
		if(it->control == -1) {
			//single qubit gate
			QASMparser::gate g = *it;
			unsigned int layer = last_layer[g.target] + 1;

			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<QASMparser::gate>());
			}
			mapped_circuit[layer].push_back(g);
			last_layer[g.target] = layer;
		} else {
			QASMparser::gate g = *it;
			unsigned int layer = max(last_layer[g.control], last_layer[g.target]) + 1;
			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<QASMparser::gate>());
			}
			mapped_circuit[layer].push_back(g);

			last_layer[g.target]  = layer;
			last_layer[g.control] = layer;
		}
	}
	
	delete[] last_layer;
}
