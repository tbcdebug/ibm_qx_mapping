#include "mapper.hpp"

#include <boost/program_options.hpp>


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

int main(int argc, char** argv) {
	bool   verbose = false;
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
			("verbose,v",       po::bool_switch(&verbose),               "verbose");

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
		ofstream of(argv[2]);

		of << "OPENQASM 2.0;" << endl;
		of << "include \"qelib1.inc\";" << endl;
		of << "qreg q[16];" << endl;
		of << "creg c[16];" << endl;

		for (vector<vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
				it != mapped_circuit.end(); it++) {
			vector<QASMparser::gate> v = *it;
			for (vector<QASMparser::gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
				of << it2->type << " ";
				if (it2->control != -1) {
					of << "q[" << it2->control << "],";
				}
				of << "q[" << it2->target << "];" << endl;
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
	nodes.delete_queue();
	
    //nodes = unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less>();
	
	node n = create_node();
	n = create_node();
	
	delete_node(n);
	n = create_node();
	delete_node(n);
	*/
	return 0;
}