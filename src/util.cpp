#include "mapper.hpp"

#include <iostream>
#include <algorithm>
#include <string.h>
#include <set>
#include <climits>
#include <fstream>

void expand_node(const vector<int>& qubits, unsigned int qubit, edge *swaps, int nswaps,
				 int* used, node base_node, const vector<QASMparser::gate>& gates, double** dist, int next_layer) {

	if (qubit == qubits.size()) {
		//base case: insert node into queue
		if (nswaps == 0) {
			return;
		}
		node new_node;

		new_node.qubits = new int[positions];
		new_node.locations = new int[nqubits];

		memcpy(new_node.qubits, base_node.qubits, sizeof(int) * positions);
		memcpy(new_node.locations, base_node.locations, sizeof(int) * nqubits);

		new_node.swaps = vector<vector<edge> >();
		new_node.nswaps = base_node.nswaps + nswaps;
		for (vector<vector<edge> >::iterator it2 = base_node.swaps.begin();
			 it2 != base_node.swaps.end(); it2++) {
            vector<edge> new_v(*it2);
			new_node.swaps.push_back(new_v);
		}

		//new_node.depth = base_node.depth + 5;
		new_node.cost_fixed = base_node.cost_fixed + 7 * nswaps;
		new_node.cost_heur = 0;

        vector<edge> new_swaps;
		for (int i = 0; i < nswaps; i++) {
			new_swaps.push_back(swaps[i]);
			int tmp_qubit1 = new_node.qubits[swaps[i].v1];
			int tmp_qubit2 = new_node.qubits[swaps[i].v2];

			new_node.qubits[swaps[i].v1] = tmp_qubit2;
			new_node.qubits[swaps[i].v2] = tmp_qubit1;

			if (tmp_qubit1 != -1) {
				new_node.locations[tmp_qubit1] = swaps[i].v2;
			}
			if (tmp_qubit2 != -1) {
				new_node.locations[tmp_qubit2] = swaps[i].v1;
			}
		}
		new_node.swaps.push_back(new_swaps);
		new_node.done = 1;

		for (vector<QASMparser::gate>::const_iterator it = gates.begin(); it != gates.end();
			 it++) {
			QASMparser::gate g = *it;
			if (g.control != -1) {
#if HEUR_ADMISSIBLE
				new_node.cost_heur = max(new_node.cost_heur, dist[new_node.locations[g.control]][new_node.locations[g.target]]);
#else
				new_node.cost_heur = new_node.cost_heur + dist[new_node.locations[g.control]][new_node.locations[g.target]];
#endif
				if(dist[new_node.locations[g.control]][new_node.locations[g.target]] > 4) {
					new_node.done = 0;
				}
			}
		}

		//Calculate heuristics for the cost of the following layer
		new_node.lookahead_penalty = 0;
#if LOOK_AHEAD
		if(next_layer != -1) {
			for (vector<QASMparser::gate>::const_iterator it = layers[next_layer].begin(); it != layers[next_layer].end();
							it++) {
                QASMparser::gate g = *it;
				if (g.control != -1) {
					if(new_node.locations[g.control] == -1 && new_node.locations[g.target] == -1) {
						//No additional penalty in heuristics
					} else if(new_node.locations[g.control] == -1) {
						int min = 1000;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[i][new_node.locations[g.target]] < min) {
								min = dist[i][new_node.locations[g.target]];
							}
						}
						new_node.lookahead_penalty = new_node.lookahead_penalty + min;
					} else if(new_node.locations[g.target] == -1) {
						int min = 1000;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[new_node.locations[g.control]][i] < min) {
								min = dist[new_node.locations[g.control]][i];
							}
						}
						new_node.lookahead_penalty = new_node.lookahead_penalty + min;
					} else {
#if HEURISTIC_ADMISSIBLE
						new_node.lookahead_penalty = max(new_node.lookahead_penalty, dist[new_node.locations[g.control]][new_node.locations[g.target]]);
#else
						new_node.lookahead_penalty = new_node.lookahead_penalty + dist[new_node.locations[g.control]][new_node.locations[g.target]];
#endif
					}
				}
			}
		}
#endif

		nodes.push(new_node);
	} else {
		expand_node(qubits, qubit + 1, swaps, nswaps, used, base_node, gates,
					dist, next_layer);

		for (set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
			edge e = *it;
			if (e.v1 == base_node.locations[qubits[qubit]]
				|| e.v2 == base_node.locations[qubits[qubit]]) {
				if (!used[e.v1] && !used[e.v2]) {
					used[e.v1] = 1;
					used[e.v2] = 1;
					swaps[nswaps].v1 = e.v1;
					swaps[nswaps].v2 = e.v2;
					expand_node(qubits, qubit + 1, swaps, nswaps + 1, used,
								base_node, gates, dist, next_layer);
					used[e.v1] = 0;
					used[e.v2] = 0;
				}
			}
		}
	}
}

node a_star_fixlayer(int layer, int* map, int* loc, double** dist) {

	int next_layer = get_next_layer(layer);

	node n;
	n.cost_fixed = 0;
	n.cost_heur = n.lookahead_penalty = 0;
	n.qubits = new int[positions];
	n.locations = new int[nqubits];
	n.swaps = vector<vector<edge> >();
	n.done = 1;

    vector<QASMparser::gate> v = vector<QASMparser::gate>(layers[layer]);
    vector<int> considered_qubits;

	//Find a mapping for all logical qubits in the CNOTs of the layer that are not yet mapped
	for (vector<QASMparser::gate>::iterator it = v.begin(); it != v.end(); it++) {
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
				int min = 1000;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[i][loc[g.target]] < min) {
						min = dist[i][loc[g.target]];
						min_pos = i;
					}
				}
				map[min_pos] = g.control;
				loc[g.control] = min_pos;
			} else if(loc[g.target] == -1) {
				int min = 1000;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[loc[g.control]][i] < min) {
						min = dist[loc[g.control]][i];
						min_pos = i;
					}
				}
				map[min_pos] = g.target;
				loc[g.target] = min_pos;
			}
			n.cost_heur = max(n.cost_heur, dist[loc[g.control]][loc[g.target]]);
		} else {
			//Nothing to do here
		}
	}

	if(n.cost_heur > 4) {
		n.done = 0;
	}

    memcpy(n.qubits, map, sizeof(int) * positions);
    memcpy(n.locations, loc, sizeof(int) * nqubits);

	nodes.push(n);

	int *used = new int[positions];
	for (int i = 0; i < positions; i++) {
		used[i] = 0;
	}
	edge *edges = new edge[considered_qubits.size()];

	//Perform an A* search to find the cheapest permutation
	while (!nodes.top().done) {
		node n = nodes.top();
		nodes.pop();

		expand_node(considered_qubits, 0, edges, 0, used, n, v, dist, next_layer);

		delete[] n.locations;
		delete[] n.qubits;
	}

	node result = nodes.top();
	nodes.pop();

	//clean up
	delete[] used;
	delete[] edges;
	while (!nodes.empty()) {
		node n = nodes.top();
		nodes.pop();
		delete[] n.locations;
		delete[] n.qubits;
	}
	return result;
}

// fix the position of the single qubit gates
void fix_positions_of_single_qubit_gates(int* locations, int* qubits, vector<QASMparser::gate>& all_gates) {
	for(vector<QASMparser::gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];
			qubits[it->control] = tmp_qubit2;
			qubits[it->target] = tmp_qubit1;

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
	for(int i=0; i<positions; i++) {
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

			last_layer[g.target] = layer;
			last_layer[g.control] = layer;
		}
	}
	
	delete[] last_layer;
}

void initial_mapping(circuit_properties& properties) {
	int* qubits    = properties.qubits;
	int* locations = properties.locations;

	for (vector<QASMparser::gate>::iterator it = layers[0].begin(); it != layers[0].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			for(set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
				if(qubits[it->v1] == -1 && qubits[it->v2] == -1) {
					qubits[it->v1] = g.control;
					qubits[it->v2] = g.target;
					locations[g.control] = it->v1;
					locations[g.target] = it->v2;
					break;
				}
			}
		}
	}
	for(unsigned int i=0; i<nqubits; i++) {
		if(locations[i] == -1) {
			int j=0;
			while(qubits[j]!=-1){
				j++;
			}
			locations[i] = j;
			qubits[j] = i;
		}
	}
}

void mapping(const vector<QASMparser::gate>& gates, vector<vector<QASMparser::gate>>& mapped_circuit, 
			vector<QASMparser::gate>& all_gates, int &total_swaps, circuit_properties& properties) {
	int* qubits    = properties.qubits;
	int* locations = properties.locations;
	
	// init layers
	layers = init_layers(gates);

#if USE_INITIAL_MAPPING
	initial_mapping(properties);	
#endif

	//Fix the mapping of each layer
	for (unsigned int i = 0; i < layers.size(); i++) {
		node result = a_star_fixlayer(i, qubits, locations, dist);

		adapt_circuit_properties(properties, result);	
	    qubits    = properties.qubits;
	    locations = properties.locations;


        vector<QASMparser::gate> h_gates = vector<QASMparser::gate>();
		//The first layer does not require a permutation of the qubits
		if (i != 0) {
			//Add the required SWAPs to the circuits
			for (vector<vector<edge> >::iterator it = result.swaps.begin();
				 it != result.swaps.end(); it++) {
				for (vector<edge>::iterator it2 = it->begin(); it2 != it->end(); it2++) {

					edge e = *it2;
					QASMparser::gate cnot;
					QASMparser::gate h1;
					QASMparser::gate h2;
					if (graph.find(e) != graph.end()) {
						cnot.control = e.v1;
						cnot.target = e.v2;
					} else {
						cnot.control = e.v2;
						cnot.target = e.v1;

						int tmp = e.v1;
						e.v1 = e.v2;
						e.v2 = tmp;
						if (graph.find(e) == graph.end()) {
                            cerr << "ERROR: invalid SWAP gate" << endl;
                            exit(2);
						}
					}
                    strcpy(cnot.type, "CX");
                    strcpy(h1.type, "U(pi/2,0,pi)");
                    strcpy(h2.type, "U(pi/2,0,pi)");
					h1.control = h2.control = -1;
					h1.target = e.v1;
					h2.target = e.v2;

					QASMparser::gate gg;
					gg.control = cnot.control;
					gg.target = cnot.target;
					strcpy(gg.type, "SWP");

					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					//Insert a dummy SWAP gate to allow for tracking the positions of the logical qubits
					all_gates.push_back(gg);
					total_swaps++;
				}
			}
		}
		
		//Add all gates of the layer to the circuit
        vector<QASMparser::gate> layer_vec = layers[i];
		for (vector<QASMparser::gate>::iterator it = layer_vec.begin();
			 it != layer_vec.end(); it++) {
			QASMparser::gate g = *it;
			if (g.control == -1) {
				//single qubit gate
				if(locations[g.target] == -1) {
					//handle the case that the qubit is not yet mapped. This happens if the qubit has not yet occurred in a CNOT gate
					QASMparser::gate g2 = g;
					g2.target = -g.target -1;
					all_gates.push_back(g2);
				} else {
					//Add the gate to the circuit
					g.target = locations[g.target];
					all_gates.push_back(g);
				}
			} else {
				//CNOT gate
				g.target  = locations[g.target];
				g.control = locations[g.control];

				edge e;
				e.v1 = g.control;
				e.v2 = g.target;

				if (graph.find(e) == graph.end()) {
					//flip the direction of the CNOT by inserting H gates
					e.v1 = g.target;
					e.v2 = g.control;
					if (graph.find(e) == graph.end()) {
                        cerr << "ERROR: invalid CNOT: " << e.v1 << " - " << e.v2 << endl;
                        exit(3);
					}
					QASMparser::gate h;
					h.control = -1;
					strcpy(h.type, "U(pi/2,0,pi)");
					h.target = g.target;
					all_gates.push_back(h);

					h_gates.push_back(h);
					h.target = g.control;
					all_gates.push_back(h);

					h_gates.push_back(h);
					int tmp = g.target;
					g.target = g.control;
					g.control = tmp;
				}
				all_gates.push_back(g);
			}
		}
		
		if (h_gates.size() != 0) {
			if (result.cost_heur == 0) {
                cerr << "ERROR: invalid heuristic cost!" << endl;
                exit(2);
			}

			for (vector<QASMparser::gate>::iterator it = h_gates.begin();
				 it != h_gates.end(); it++) {
				all_gates.push_back(*it);
			}
		}

	}	

	fix_positions_of_single_qubit_gates(locations, qubits, all_gates);
	generate_circuit(mapped_circuit, all_gates);
}
