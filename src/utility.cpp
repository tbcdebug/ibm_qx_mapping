#include "mapper.hpp"

// checks if a vector contains a certain edge
bool contains(const vector<int>& v, const int e) {
	for (vector<int>::const_iterator it = v.begin(); it != v.end(); it++) {
		if (*it == e) {
			return true;
		}
	}
	return false;
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

// generates a circuit based on the gates
void generate_circuit(vector<vector<QASMparser::gate>>& mapped_circuit, const vector<QASMparser::gate>& all_gates) {
    int *last_layer = new int[positions];
	for(int i=0; i<positions; i++) {
		last_layer[i] = -1;
	}

	//build resulting circuit
	for(vector<QASMparser::gate>::const_iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
#if VERIFICATION
			QASMparser::gate g;
			g.control = it->control;
			g.target  = it->target;
			strcpy(g.type, "SWAP");
			unsigned int layer = max(last_layer[g.control], last_layer[g.target]) + 1;
			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<QASMparser::gate>());
			}
			mapped_circuit[layer].push_back(g);
			last_layer[g.target]  = layer;
			last_layer[g.control] = layer;
#endif
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

//Fix the position of the single qubit gates 
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
				int loc;
				for(loc = 0; qubits[loc] != -1; loc++); //DONE for loop
				locations[target] = loc;
			}
		}
	}
}