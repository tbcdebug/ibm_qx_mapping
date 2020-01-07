<<<<<<< HEAD
#include <iostream>
#include <algorithm>
#include <string.h>
#include <set>
#include <climits>
#include <fstream>
#include <regex>
#include <math.h>

#include "QASMparser.h"
#include "unique_priority_queue.h"

#define LOOK_AHEAD 1
#define HEURISTIC_ADMISSIBLE 0
#define USE_INITIAL_MAPPING 0
#define MINIMAL_OUTPUT 1
#define DUMP_MAPPED_CIRCUIT 1
#define DUMP_STATISTICS 1

#define ARCH_LINEAR_N 0
#define ARCH_IBM_QX5 1

#ifndef ARCH
// assume default architecture
#define ARCH ARCH_LINEAR_N
#endif

using namespace std;

int** dist;
int positions;
unsigned long ngates = 0;
unsigned int nqubits = 0;

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

struct node {
	int cost_fixed;
	int cost_heur;
	int cost_heur2;
	int depth;
	int* qubits; // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int* locations; // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
	int nswaps;
	int done;
	std::vector<std::vector<edge>> swaps;
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
		if ((x.cost_fixed + x.cost_heur + x.cost_heur2) != (y.cost_fixed + y.cost_heur + y.cost_heur2)) {
			return (x.cost_fixed + x.cost_heur + x.cost_heur2) > (y.cost_fixed + y.cost_heur + y.cost_heur2);
		}
=======
#include "mapper.hpp"
>>>>>>> 3b65bb658d3b0ce0953fa64f3f095a245c407093

#include <boost/program_options.hpp>
#include <regex>
#include <math.h>

namespace po = boost::program_options;

/**
 * Global variables
 */
double**       dist;
int            positions;
unsigned long  ngates  = 0;
unsigned int   nqubits = 0;

set<edge>                                                                    graph;
vector<vector<QASMparser::gate>>                                             layers;
unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less> nodes;

double get_pi_div(double val) {
	if(val == 0) {
		return 0;
	}
	const int precision = 10000;
	return round(M_PI / val * precision) / precision;
}

int main(int argc, char** argv) {
	bool   verbose     = false;
	bool   real_format = false;
	string input,  input_coupling;
	string output, output_statistics;

	// argument handling
	try {
		po::options_description desc{"Options"};
    	desc.add_options()
			("help,h",                                                   "help screen")
			("input,i",         po::value<string>(&input)->required(),   "input file")
			("output,o",        po::value<string>(&output),              "output file")
			("statistic,s",     po::value<string>(&output_statistics),   "output statistics file")
			("coupling_file,c", po::value<string>(&input_coupling),      "coupling graph - file")
			("verbose,v",       po::bool_switch(&verbose),               "verbose")
			("real,r",          po::bool_switch(&real_format),           "output the circuit in the real format");

		po::positional_options_description p;
		p.add("input",     1);
		p.add("statistic", 1);
		p.add("output",    1);

		po::variables_map vm; 
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		
		if (vm.count("help")) {
        	cout << desc << endl;
            return 0;
        }
		po::notify(vm);
	} catch (const po::error &ex) {
		cerr << ex.what() << endl;
		exit(ERROR);
	}

	if(verbose) {
		cout << "Input:        " << input             << endl;
		cout << "Output:       " << output            << endl;
		cout << "Statistic:    " << output_statistics << endl;
		cout << "CouplingFile: " << input_coupling    << endl;
		cout << "Verbose:      " << verbose           << endl;
	}

	// parsing	
	QASMparser* parser = new QASMparser(input.c_str());
	parser->Parse();

	vector<QASMparser::gate> gates = parser->getGates();
	nqubits = parser->getNqubits();
	ngates  = parser->getNgates();

	parser->clear();
	delete parser;

	
	// graph handling
	if(!generate_graph(input_coupling)) {
		cout << "Error while generating the graph" << endl;
		exit(ERROR);
	}
	if(positions > 0 ? nqubits > (unsigned int)positions : (int)nqubits > positions) {
        cerr << "ERROR before mapping: more logical qubits than physical ones!" << endl;
        exit(ERROR);
    }

		
	// print infos
	const char* bName = basename(input.c_str());
	if(verbose) {
		cout << "Circuit name: " << bName << " (requires " << nqubits << " qubits)" << endl;
		cout << endl;
		cout << "Before mapping: "                      << endl;
		cout << "  elementary gates: " << ngates        << endl;
		cout << "  depth:            " << layers.size() << endl;
	} else {
    	cout << bName << ',' << nqubits << ',' << ngates << ',' << layers.size() << ',' << flush;
	}
	
	// start mapping algorithm
	clock_t begin_time = clock();

	int                              total_swaps = 0;	
	circuit_properties               properties  = create_circuit_properties();
    vector<QASMparser::gate>         all_gates;
	vector<vector<QASMparser::gate>> mapped_circuit;

	mapping(gates, mapped_circuit, all_gates, total_swaps, properties);

	double    time     = double(clock() - begin_time) / CLOCKS_PER_SEC;
	int       depth    = mapped_circuit.size();
	int       cost     = all_gates.size()-total_swaps;
#if SPECIAL_OPT	
	long long fidelity = fidelity_cost(properties.fidelities);
#else
	long long fidelity = 0;
#endif
	// print statistics
	if(verbose) {
		cout << endl << "After mapping (no post mapping optimizations are conducted): " << endl;
		cout << "  elementary gates: " << cost  << endl;
		cout << "  depth:            " << depth << endl;

		cout << "\nThe mapping required " << time << " seconds" << endl;

		cout << "\nInitial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << endl;

		for(uint i = 0; i < nqubits; i++) {
			cout << "  q" << i << " is initially mapped to Q" << properties.locations[i] << endl;
		} 
	} else {
    	cout << time << ',' << cost << ',' << depth << endl;
	}

	// dump resulting circuit
	if(!output.empty()) {
		ofstream of(output);
		if(real_format) {
			of << ".numvars "   << nqubits << endl;
			of << ".variables";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << " q" << i;
			}
			of << endl;
			of << ".constants ";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << "0";
			}
			of << endl;
			of << ".begin" << endl;
			for (vector<vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
					it != mapped_circuit.end(); it++) {
				vector<QASMparser::gate> v = *it;
				for (vector<QASMparser::gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
					string hadamard = "U(pi/2,0,pi)";
					if(it2->control != -1) {
						of << "t2 " << "q" << it2->control << " q" << it2->target << endl;
					} else if(hadamard.compare(it2->type) == 0) {
						of << "h1 q" << it2->target << endl;
					} else {
						std::string s(it2->type);
						std::regex rgx("U\\(([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+)\\).*");
						std::smatch match;
						
						if(std::regex_search(s, match, rgx)) {							
							double theta = stof(match[1]);
							double phi   = stof(match[3]);
							double delta = stof(match[5]);

							double theta_div = get_pi_div(theta); 
							double phi_div   = get_pi_div(phi); 
							double delta_div = get_pi_div(delta);

							/*
							cout << "THEATA" << theta << endl;
							cout << "DIV   " << theta_div << endl;
							cout << "PHI   " << phi << endl;
							cout << "DIV   " << phi_div << endl;
							cout << "DELTA " << delta << endl;
							cout << "DIV   " << delta_div << endl;
							*/

							// conversion to rotation gates
							if(phi_div == 0) {
								of << "rz1:" << 1                     << " q" << it2->target << endl; //1.0 / 3
							} else {
								of << "rz1:" << (int)(phi_div / (1 + 3 * phi_div)) << " q" << it2->target << endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << endl;
							if(theta_div == 0) {
								of << "rz1:" << 1                           << " q" << it2->target << endl;
							} else {
								of << "rz1:" << (int)(theta_div / (1 + theta_div)) << " q" << it2->target << endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << endl;
							if(delta_div != 0) {
								of << "rz1:" << delta_div                   << " q" << it2->target << endl;
							}
							//for(auto x: match)
							//	std::cout << "match: " << x << '\n';
							//s = match.suffix().str();
						}
					}
				}
			}
			of << ".end" << endl;
		} else {
			of << "OPENQASM 2.0;"              << endl;
			of << "include \"qelib1.inc\";"    << endl;
			of << "qreg q[" << nqubits << "];" << endl;
			of << "creg c[" << nqubits << "];" << endl;

			for (vector<vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
				it != mapped_circuit.end(); it++) {
				for (vector<QASMparser::gate>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
					of << it2->type << " ";
					if (it2->control != -1) {
						of << "q[" << it2->control << "],";
					}
					of << "q[" << it2->target << "];" << endl;
				}
			}
		}
	}

	// store timing
	if(!output_statistics.empty()) {
		ofstream ofstat (output_statistics, ofstream::app);
		//ofstat << bName << " : " << time << " " << depth << " " << cost << " " << fidelity << " " << alloc_tries << " " << total_swaps << endl;
		ofstat << bName << " : " << time << " " << depth << " " << cost << " " << fidelity << " " << total_swaps << endl;
	}

	delete_circuit_properties(properties);
	
	
	/*
	positions = 16;
	nqubits   = 16;
	
	for(int i = 0; i < 2; i++) {
		node n = create_node();
		//Initially, no physical qubit is occupied
		for (int j = 0; j < positions; j++) {
			n.qubits[j] = -1;
		}

		//Initially, no logical qubit is mapped to a physical one
		for(unsigned j = 0; j < nqubits; j++) {
			n.locations[j] = -1;
		}
		nodes.push(n);
	}
<<<<<<< HEAD
	return result;
}

double get_pi_div(double val) {
	if(val == 0) {
		return 0;
	}
	const int precision = 10000;
	return round(M_PI / val * precision) / precision;
}

int main(int argc, char** argv) {
/*
#if DUMP_MAPPED_CIRCUIT
	if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        std::exit(1);
=======
	nodes.delete_queue();
	
    //nodes = unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less>();
	
	node n = create_node();
	n = create_node();
	
	delete_node(n);
	n = create_node();
	delete_node(n);
	*/

	/*
	std::string s = "u(0.90, 0.00, -0.7)";
    std::regex rgx("u\\(([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+)\\).*");
    
	std::smatch match;
	
    if(std::regex_search(s, match, rgx)) {
        double theta = stof(match[1]);
		double phi   = stof(match[3]);
		double delta = stof(match[5]);
		
		cout << theta << endl;
		cout << phi   << endl;
		cout << delta << endl;
		//for(auto x: match)
		//	std::cout << "match: " << x << '\n';
		//s = match.suffix().str();
>>>>>>> 3b65bb658d3b0ce0953fa64f3f095a245c407093
	}
	*/
	/*bool   verbose     = false;
	bool   real_format = false;
	string input,  input_coupling;
	string output, output_statistics;

	// argument handling
	try {
		po::options_description desc{"Options"};
    	desc.add_options()
			("help,h",                                                   "help screen")
			("input,i",         po::value<string>(&input)->required(),   "input file")
			("output,o",        po::value<string>(&output),              "output file")
			("statistic,s",     po::value<string>(&output_statistics),   "output statistics file")
			("coupling_file,c", po::value<string>(&input_coupling),      "coupling graph - file")
			("verbose,v",       po::bool_switch(&verbose),               "verbose")
			("real,r",          po::bool_switch(&real_format),           "output the circuit in the real format");

		po::positional_options_description p;
		p.add("input",     1);
		p.add("statistic", 1);

		po::variables_map vm; 
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		
		if (vm.count("help")) {
        	cout << desc << endl;
            return 0;
        }
		po::notify(vm);
	} catch (const po::error &ex) {
		cerr << ex.what() << endl;
		exit(ERROR);
	}
<<<<<<< HEAD
#endif
*/
	QASMparser* parser = new QASMparser(argv[1]);
=======

	// parsing	
	QASMparser* parser = new QASMparser(input.c_str());
>>>>>>> 3b65bb658d3b0ce0953fa64f3f095a245c407093
	parser->Parse();

	vector<QASMparser::gate> gates = parser->getGates();
	nqubits = parser->getNqubits();
	ngates  = parser->getNgates();

	parser->clear();
	delete parser;
	real_format = true;

	// dump resulting circuit
	if(!output.empty()) {
		ofstream of(output);
		if(real_format) {
			of << ".numvars "   << nqubits << endl;
			of << ".variables";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << " q" << i;
			}
			of << endl;
			of << ".constants ";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << "0";
			}
			of << endl;
			of << ".begin" << endl;
			for (vector<QASMparser::gate>::iterator it = gates.begin();
					it != gates.end(); it++) {
				string hadamard = "U(pi/2,0,pi)";
				if(it->control != -1) {
					of << "t2 " << "q[" << it->control << "] q[" << it->target << "];" << endl;
				} else if(hadamard.compare(it->type) == 0) {
					of << "h1 q[" << it->target << "];" << endl;
				} else {
					std::string s(it->type);
					std::regex rgx("u\\(([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+)\\).*");
					
					std::smatch match;
					if(std::regex_search(s, match, rgx)) {
						double theta = stof(match[1]);
						double phi   = stof(match[3]);
						double delta = stof(match[5]);
						
						cout << theta << endl;
						cout << phi   << endl;
						cout << delta << endl;
						//for(auto x: match)
						//	std::cout << "match: " << x << '\n';
						//s = match.suffix().str();
					}
				}
			}
			of << ".end" << endl;
		} else {
		}
	}
<<<<<<< HEAD

	double time = double(clock() - begin_time) / CLOCKS_PER_SEC;

#if !MINIMAL_OUTPUT
    std::cout << std::endl << "After mapping (no post mapping optimizations are conducted): " << std::endl;
	std::cout << "  elementary gates: " << all_gates.size()-total_swaps << std::endl;
	std::cout << "  depth: " << mapped_circuit.size() << std::endl;

	std::cout << "\nThe mapping required " << time << " seconds" << std::endl;

	std::cout << "\nInitial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << std::endl;

	for(int i=0; i<nqubits; i++) {
		std::cout << "  q" << i << " is initially mapped to Q" << locations[i] << std::endl;
	}
#else
    std::cout << time << ',' << (all_gates.size()-total_swaps) << ',' << mapped_circuit.size() << std::endl;
#endif

#if DUMP_MAPPED_CIRCUIT
	assert(argc > 3);
	//Dump resulting circuit
	std::ofstream of(argv[3]);
	bool real_format = false;
	if(real_format) {
			of << ".numvars "   << nqubits << endl;
			of << ".variables";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << " q" << i;
			}
			of << endl;
			of << ".constants ";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << "0";
			}
			of << endl;
			of << ".begin" << endl;
			for (vector<vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
					it != mapped_circuit.end(); it++) {
				vector<QASMparser::gate> v = *it;
				for (vector<QASMparser::gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
					string hadamard = "U(pi/2,0,pi)";
					if(it2->control != -1) {
						of << "t2 " << "q" << it2->control << " q" << it2->target << endl;
					} else if(hadamard.compare(it2->type) == 0) {
						of << "h1 q" << it2->target << endl;
					} else {
						std::string s(it2->type);
						std::regex rgx("U\\(([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+)\\).*");
						std::smatch match;
						
						if(std::regex_search(s, match, rgx)) {							
							double theta = stof(match[1]);
							double phi   = stof(match[3]);
							double delta = stof(match[5]);

							double theta_div = get_pi_div(theta); 
							double phi_div   = get_pi_div(phi); 
							double delta_div = get_pi_div(delta);

							/*
							cout << "THEATA" << theta << endl;
							cout << "DIV   " << theta_div << endl;
							cout << "PHI   " << phi << endl;
							cout << "DIV   " << phi_div << endl;
							cout << "DELTA " << delta << endl;
							cout << "DIV   " << delta_div << endl;
							*/

							// conversion to rotation gates
							if(phi_div == 0) {
								of << "rz1:" << 1                     << " q" << it2->target << endl; //1.0 / 3
							} else {
								of << "rz1:" << (int)(phi_div / (1 + 3 * phi_div)) << " q" << it2->target << endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << endl;
							if(theta_div == 0) {
								of << "rz1:" << 1                           << " q" << it2->target << endl;
							} else {
								of << "rz1:" << (int)(theta_div / (1 + theta_div)) << " q" << it2->target << endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << endl;
							if(delta_div != 0) {
								of << "rz1:" << delta_div                   << " q" << it2->target << endl;
							}
							//for(auto x: match)
							//	std::cout << "match: " << x << '\n';
							//s = match.suffix().str();
						}
					}
				}
			}
			of << ".end" << endl;
		} else {
			of << "OPENQASM 2.0;" << endl;
			of << "include \"qelib1.inc\";" << endl;
			of << "qreg q[16];" << endl;
			of << "creg c[16];" << endl;

			for (vector<vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
				it != mapped_circuit.end(); it++) {
				for (vector<QASMparser::gate>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
					of << it2->type << " ";
					if (it2->control != -1) {
						of << "q[" << it2->control << "],";
					}
					of << "q[" << it2->target << "];" << endl;
				}
			}
		}
#endif
#if DUMP_STATISTICS
    assert(argc > 2);
	int       depth    = mapped_circuit.size();
	int       cost     = all_gates.size()-total_swaps;
	std::ofstream ofstat (argv[2], std::ofstream::app);
	//ofstat << bName << " : " << time << " " << depth << " " << cost << " " << fidelity << " " << alloc_tries << " " << total_swaps << endl;
	ofstat << bName << " : " << time << " " << depth << " " << cost << " " << 0 << " " << " " << total_swaps << std::endl;
#endif
	delete[] locations;
	delete[] qubits;
	delete[] last_layer;

=======
	*/
>>>>>>> 3b65bb658d3b0ce0953fa64f3f095a245c407093
	return 0;
}