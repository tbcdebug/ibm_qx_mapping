#include "mapper.hpp"

#include <boost/program_options.hpp>


namespace po = boost::program_options;

/**
 * Global variables
 */
double**       dist;
int            positions;
unsigned long  ngates = 0;
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
			("help,h",                               "help screen")
			("input,i",         po::value<string>(), "input file")
			("output,o",        po::value<string>(), "output file")
			("statistic,s",     po::value<string>(), "output statistics file")
			("coupling_file,c", po::value<string>(), "coupling graph - file")
			("verbose,v",                            "verbose");

		po::positional_options_description p;
		p.add("input", -1);

		po::variables_map vm; 
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);

		if(vm.count("help")) {
			cout << desc << endl;
			exit(SUCCESS);
		}
		if(vm.count("input")) {
			input = vm["input"].as<string>();
		} else {			
			cerr << "Input has to be specified." << endl;
			exit(ERROR);
		}
		if(vm.count("output")) {
			output = vm["output"].as<string>();
		}
		if(vm.count("statistic")) {
			output_statistics = vm["statistic"].as<string>();
		}
		if(vm.count("coupling_file")) {
			input_coupling = vm["coupling_file"].as<string>();
		}
		if(vm.count("verbose")) {
			verbose = true;
		}

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

		cout << endl << "Before mapping: " << endl;
		cout << "  elementary gates: " << ngates << endl;
		cout << "  depth: " << layers.size() << endl;
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

	double time = double(clock() - begin_time) / CLOCKS_PER_SEC;
	
	// print statistics
	if(verbose) {
		cout << endl << "After mapping (no post mapping optimizations are conducted): " << endl;
		cout << "  elementary gates: " << all_gates.size()-total_swaps << endl;
		cout << "  depth: " << mapped_circuit.size() << endl;

		cout << "\nThe mapping required " << time << " seconds" << endl;

		cout << "\nInitial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << endl;

		for(uint i = 0; i < nqubits; i++) {
			cout << "  q" << i << " is initially mapped to Q" << properties.locations[i] << endl;
		} 
	} else {
    	cout << time << ',' << (all_gates.size()-total_swaps) << ',' << mapped_circuit.size() << endl;
	}

	if(!output.empty()) {
		//Dump resulting circuit
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

	delete_circuit_properties(properties);
}