#include "mapper.hpp"

#include <string.h>

/**
 * static function declarations
 */
static void lookahead(const int next_layer, node& new_node);
static node a_star_fixlayer(const int layer, circuit_properties& properties);
#if ONE_SWAP_PER_EXPAND
static void expand_node_add_one_swap(const vector<int>& qubits, const edge e, node base_node, const vector<QASMparser::gate>& gates,
							  const int next_layer);
static void expand_node(const vector<int>& qubits, const node base_node, const vector<QASMparser::gate>& gates, const int next_layer);
#else
static void expand_node(const vector<int>& qubits, unsigned int qubit, edge *swaps, int nswaps,
				 int* used, node base_node, const vector<QASMparser::gate>& gates, int next_layer);
#endif

/**
 * calculates and sets the lookahead cost for the specified node
 * may consider many "future" layers based on N_LOOK_AHEADS
 */
static void lookahead(const int layer, node& new_node) {
	int    next_layer = layer;
	double factor     = FIRST_LOOK_AHEAD_FACTOR;
	for(int i = 0; i < N_LOOK_AHEADS; i++) {
		if(next_layer == -1) {
			break;
		}
		double penalty = 0;
#if SPECIAL_OPT
		// considers the depth and fidelity of the layers between this and the next layer with cnots
		if(i == 0 && (DEPTH_PERCENTAGE > 0 || FIDELITY_NORM > 0)) { // calculate depth for the first next layer
			int* depths     = new int[positions];
			int* fidelities = new int[positions];
			
			memcpy(depths,     new_node.depths,     sizeof(int) * positions);
			memcpy(fidelities, new_node.fidelities, sizeof(int) * positions);

			for(int j = layer + 1; j <= next_layer; j++) {
				for (vector<QASMparser::gate>::const_iterator it = layers[j].begin(); it != layers[j].end();
						it++) {
					QASMparser::gate g = *it;
					int target = new_node.locations[g.target];
					if (g.control != -1) {
						int control = new_node.locations[g.control];
						if(control == -1 && target == -1) {
							//No additional penalty in heuristics
						} else if(control == -1) {
							depths[target]     += DEPTH_GATE;
							fidelities[target] += FIDELITY_GATE;
						} else if(target == -1 ) {
							depths[control]     += DEPTH_GATE;
							fidelities[control] += FIDELITY_GATE;
						} else {
							if(dist[control][target] < 1) {
								depths[control]     += DEPTH_GATE;	
								depths[target]      += DEPTH_GATE;	
								fidelities[control] += FIDELITY_CNOT;	
								fidelities[target]  += FIDELITY_CNOT;		
							} else {
								depths[control]     += DEPTH_SWAP;	
								depths[target]      += DEPTH_SWAP;	
								fidelities[control] += FIDELITY_SWAP;	
								fidelities[target]  += FIDELITY_SWAP;										
							}
						}
					} else {
						depths[target]     += DEPTH_GATE;
						fidelities[target] += FIDELITY_GATE;	
					}	
				}
			}
			new_node.lookahead_penalty = factor * 
					  (get_maximal_depth(depths) / ((double) DEPTH_SWAP) * DEPTH_PERCENTAGE + 
					  fidelity_cost(fidelities)  * FIDELITY_NORM);
			
			delete[] depths;
			delete[] fidelities;			
		}
#endif
		for (vector<QASMparser::gate>::const_iterator it = layers[next_layer].begin(); it != layers[next_layer].end();
						it++) {
			QASMparser::gate g = *it;
			if (g.control != -1) {
				if(new_node.locations[g.control] == -1 && new_node.locations[g.target] == -1) {
					//No additional penalty in heuristics
				} else if(new_node.locations[g.control] == -1) {
					int min = 1000;
					for(int i = 0; i < positions; i++) {
						if(new_node.qubits[i] == -1 && dist[i][new_node.locations[g.target]] < min) {
							min = dist[i][new_node.locations[g.target]];
						}
					}
					penalty = heuristic_function(penalty, min);
				} else if(new_node.locations[g.target] == -1) {
					int min = 1000;
					for(int i = 0; i < positions; i++) {
						if(new_node.qubits[i] == -1 && dist[new_node.locations[g.control]][i] < min) {
							min = dist[new_node.locations[g.control]][i];
						}
					}
					penalty = heuristic_function(penalty, min);
				} else {
					penalty = get_heuristic_cost(penalty, new_node, g);
				}
			}
			
		}
#if SPECIAL_OPT
		new_node.lookahead_penalty += factor * penalty * COST_PERCENTAGE;
#else
		new_node.lookahead_penalty += factor * penalty;
#endif
		factor                      *= GENERAL_LOOK_AHEAD_FACTOR;
		next_layer = get_next_layer(next_layer); //TODO single qubit -> approach 
	}	
}

#if ONE_SWAP_PER_EXPAND

/**
 * creates a node based on the old node and the specified swap/edge
 */
static void expand_node_add_one_swap(const vector<int>& qubits, const edge e, node base_node, 
                                     const vector<QASMparser::gate>& gates, const int next_layer) {
	// create new node based on original one
	node new_node = create_node(base_node, &e);
	
	// calculate heuristics and check if node is target
	for (vector<QASMparser::gate>::const_iterator it = gates.begin(); it != gates.end();
		 it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			new_node.cost_heur = get_heuristic_cost(new_node.cost_heur, new_node, g);
			check_if_not_done(new_node, dist[new_node.locations[g.control]][new_node.locations[g.target]]);
		}
	}
#if LOOK_AHEAD
	lookahead(next_layer, new_node);
#endif
	nodes.push(new_node);		
}

/**
 * adds all possible swaps to the current node considering the specified qubits
 */
static void expand_node(const vector<int>& qubits, const node base_node, const vector<QASMparser::gate>& gates, const int next_layer) {
	int **used_swaps;
	used_swaps = new int*[nqubits];
	for(unsigned int i = 0; i < nqubits; i++) {
		used_swaps[i] = new int[nqubits]();
	}

	// successor = one swap for each considered qubit for each edge in the coupling graph
	for(unsigned int qubit = 0; qubit < qubits.size(); qubit++) {
		for (set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
			edge e = *it;			
			if (e.v1 == base_node.locations[qubits[qubit]] 			// swap operations only according to coupling graph
		     || e.v2 == base_node.locations[qubits[qubit]]) {
				int q1 = base_node.qubits[e.v1];
				int q2 = base_node.qubits[e.v2];
				if(q2 == -1 || q1 == -1) {
					expand_node_add_one_swap(qubits, e, base_node, gates, next_layer);
				} else if (!used_swaps[q1][q2]) { 
					used_swaps[q1][q2] = 1;
					used_swaps[q2][q1] = 1;
					expand_node_add_one_swap(qubits, e, base_node, gates, next_layer);
				}
			}
		}
	}
	
	for(unsigned int i = 0; i < nqubits; i++) {
		delete[] used_swaps[i];
	}
	delete[] used_swaps;
}
#else

/**
 * adds swaps considering all permutations and the specified qubitss
 */
static void expand_node(const vector<int>& qubits, unsigned int qubit, edge *swaps, int nswaps,
				 int* used, node base_node, const vector<QASMparser::gate>& gates, int next_layer) {

	if (qubit == qubits.size()) {
		//base case: insert node into queue
		if (nswaps == 0) {
			return;
		}
		node new_node = create_node(base_node, swaps, nswaps);

		for (vector<QASMparser::gate>::const_iterator it = gates.begin(); it != gates.end();
			 it++) {
			QASMparser::gate g = *it;
			if (g.control != -1) {
				new_node.cost_heur = get_heuristic_cost(new_node.cost_heur, new_node, g);
				check_if_not_done(new_node, dist[new_node.locations[g.control]][new_node.locations[g.target]]);
			}
		}

		//Calculate heuristics for the cost of the following layer
#if LOOK_AHEAD
		lookahead(next_layer, new_node);
#endif

		nodes.push(new_node);
	} else {
		expand_node(qubits, qubit + 1, swaps, nswaps, used, base_node, gates, next_layer);

		for (set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
			edge e = *it;
			if (e.v1 == base_node.locations[qubits[qubit]]
				|| e.v2 == base_node.locations[qubits[qubit]]) {
				if (!used[e.v1] && !used[e.v2]) {
					used[e.v1] = 1;
					used[e.v2] = 1;
					swaps[nswaps].v1 = e.v1;
					swaps[nswaps].v2 = e.v2;
					expand_node(qubits, qubit + 1, swaps, nswaps + 1, used,	base_node, gates, next_layer);
					used[e.v1] = 0;
					used[e.v2] = 0;
				}
			}
		}
	}
}
#endif

/**
 * executes and a star search for the specified layer
 */
node a_star_fixlayer(const int layer, circuit_properties& properties) {
	int  next_layer = get_next_layer(layer);
	node orig       = create_node();
    
	vector<int>                    considered_qubits;
	const vector<QASMparser::gate> gates = layers[layer];
    
	map_unmapped_gates(layer, properties, orig, considered_qubits);

	update_node(orig, properties);
	check_if_not_done(orig, orig.cost_heur);
	

	nodes.push(orig);

#if USE_QUEUE_LIMIT
	unsigned int  max_node_size = MAX_QUEUE_SIZE;
#else 
	unsigned int  max_node_size = 0;
#endif

#if ONE_SWAP_PER_EXPAND
	while (!nodes.top().done) {
		node n = nodes.top();
		nodes.pop();
		
		try {
			if(max_node_size != 0 && nodes.size() > max_node_size) {
				nodes.update();
			}
			expand_node(considered_qubits, n, gates, next_layer);
		} catch (std::bad_alloc &e) {
			cout << "bad_alloc catched" << endl;
			try {
				max_node_size = nodes.size() - MAX_NODES_MARGIN;
				nodes.update();				
			} catch (std::bad_alloc &e) {
				cerr << endl;
				cerr << "bad_alloc while reseting queue -> reset to layer begin" << endl;
				nodes.restart(orig);
			}
		}

		delete_node(n);
	}
#else
	int *used = new int[positions];
	for (int i = 0; i < positions; i++) {
		used[i] = 0;
	}
	edge *edges = new edge[considered_qubits.size()];
	//Perform an A* search to find the cheapest permutation
	while (!nodes.top().done) {
		node n = nodes.top();
		nodes.pop();

		//expand_node(considered_qubits, 0, edges, 0, used, n, gates, next_layer);
		try {
			if(max_node_size != 0 && nodes.size() > max_node_size) {
				nodes.update();
			}
			expand_node(considered_qubits, 0, edges, 0, used, n, gates, next_layer);
		} catch (std::bad_alloc &e) {
			cout << "bad_alloc catched" << endl;
			try {
				max_node_size = nodes.size() - MAX_NODES_MARGIN;
				nodes.update();				
			} catch (std::bad_alloc &e) {
				cerr << endl;
				cerr << "bad_alloc while reseting queue -> reset to layer begin" << endl;
				nodes.restart(orig);
			}
		}

		delete_node(n);
	}
	//clean up
	delete[] used;
	delete[] edges;
#endif



	node result = nodes.top();
	nodes.pop();


	delete_nodes();
	return result;
}

/**
 * executes the whole mapping process
 */
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
		node result = a_star_fixlayer(i, properties);

		adapt_circuit_properties(properties, result);	
#if SPECIAL_OPT
		update_properties(properties, i);
#endif 
	    qubits    = properties.qubits;
	    locations = properties.locations;

        vector<QASMparser::gate> h_gates = vector<QASMparser::gate>();
		//The first layer does not require a permutation of the qubits
		if (i != 0) {
			//Add the required SWAPs to the circuits
			for (SWAP_LIST_TYPE::iterator it = result.swaps.begin();
				 it != result.swaps.end(); it++) {
#if ONE_SWAP_PER_EXPAND
				 {
					 edge e = *it;
#else
				for (vector<edge>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
					edge e = *it2;
#endif				
				
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
                    strcpy(h1.type,   "U(pi/2,0,pi)");
                    strcpy(h2.type,   "U(pi/2,0,pi)");
					h1.control = h2.control = -1;
					h1.target  = e.v1;
					h2.target  = e.v2;

					QASMparser::gate gg;
					gg.control = cnot.control;
					gg.target  = cnot.target;
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
					int tmp   = g.target;
					g.target  = g.control;
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
