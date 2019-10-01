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
		cout << "Output:       " << input_coupling    << endl;
		cout << "Statistic:    " << output            << endl;
		cout << "CouplingFile: " << output_statistics << endl;
		cout << "Verbose:      " << verbose           << endl;
	}

	// parsing	
	QASMparser* parser = new QASMparser(argv[1]);
	parser->Parse();

	layers  = parser->getLayers();
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

	mapper(argc, argv);
}