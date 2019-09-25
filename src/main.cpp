#include "mapper.hpp"

#include <boost/program_options.hpp>


using namespace std;
namespace po = boost::program_options;

/*
 * definition of global variables
 */
double**      dist;
int           positions;
unsigned long ngates  = 0;
unsigned int  nqubits = 0;

#if USE_QUEUE_LIMIT
unsigned int  max_node_size = MAX_QUEUE_SIZE;
#else 
unsigned int  max_node_size = 0;
#endif

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
		cout << "Output:       " << input_coupling    << endl;
		cout << "Statistic:    " << output            << endl;
		cout << "CouplingFile: " << output_statistics << endl;
		cout << "Verbose:      " << verbose           << endl;
	}

	// parse qasm
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
		std::cout << "Circuit name: " << bName << " (requires " << nqubits << " qubits)" << std::endl;

		std::cout << std::endl << "Before mapping: " << std::endl;
		std::cout << "  elementary gates: " << ngates << std::endl;
		std::cout << "  depth: " << layers.size() << std::endl;
	} else {
    	std::cout << bName << ',' << nqubits << ',' << ngates << ',' << layers.size() << ',' << std::flush;
	}

	// start mapping
	int                              total_swaps = 0;
	circuit_properties               properties  = create_circuit_properties();
	vector<vector<QASMparser::gate>> mapped_circuit;
	vector<QASMparser::gate>         all_gates;
	
	clock_t begin_time = clock();
	mapper(gates, mapped_circuit, all_gates, total_swaps, properties);
	double time = double(clock() - begin_time) / CLOCKS_PER_SEC;
	
	int *locations = properties.locations;

	if(verbose) {
		std::cout << std::endl << "After mapping (no post mapping optimizations are conducted): " << std::endl;
		std::cout << "  elementary gates: " << all_gates.size()-total_swaps << std::endl;
		std::cout << "  depth: " << mapped_circuit.size() << std::endl;

		std::cout << "\nThe mapping required " << time << " seconds" << std::endl;

		std::cout << "\nInitial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << std::endl;

		for(unsigned int i = 0; i < nqubits; i++) {
			std::cout << "  q" << i << " is initially mapped to Q" << locations[i] << std::endl;
		}
	} else {
    	std::cout << time << ',' << (all_gates.size()-total_swaps) << ',' << mapped_circuit.size() << std::endl;
	}
	
	//Dump resulting circuit
	if(!output.empty()) {
		std::ofstream of(output);

		of << "OPENQASM 2.0;"                << std::endl;
		of << "include \"qelib1.inc\";"      << std::endl;
		of << "qreg q[" << positions << "];" << std::endl;
		of << "creg c[" << positions << "];" << std::endl;

		for (std::vector<std::vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
				it != mapped_circuit.end(); it++) {
			std::vector<QASMparser::gate> v = *it;
			for (std::vector<QASMparser::gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
				of << it2->type << " ";
				if (it2->control != -1) {
					of << "q[" << it2->control << "],";
				}
				of << "q[" << it2->target << "];" << std::endl;
			}
		}
	}

	// store statistics
	if(!output_statistics.empty()) {
		ofstream ofstat(output_statistics, ofstream::app);
		//ofstat << bName << " : " << time << " " << depth << " " << cost << " " << fidelity << " " << alloc_tries << " " << total_swaps << endl; //TODO
	}

	delete_circuit_properties(properties);

	exit(SUCCESS);
}
