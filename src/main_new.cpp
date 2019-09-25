#include "mapper.hpp"


#include <iostream>
#include <boost/program_options.hpp>


using namespace std;

#define SUCCESS 0
#define ERROR   1

namespace po = boost::program_options;

int main(int argc, char** argv) {
	bool   verbose   = false;
	bool   no_output = false;
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
			("verbose,v",                            "verbose")
			("no_output,n",                          "no outsput");

		po::positional_options_description p;
		p.add("input", -1);

		po::variables_map vm; 
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);

		if(vm.count("help")) {
			cout << desc << endl;
			return SUCCESS;
		}
		if(vm.count("input")) {
			input = vm["input"].as<string>();
		} else {			
			cerr << "Input has to be specified." << endl;
			return ERROR;
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
		if(vm.count("no_output")) {
			no_output = true;
		}

	} catch (const po::error &ex) {
		cerr << ex.what() << endl;
		return ERROR;
	}

	// graph handling
	if(!generate_graph(input_coup)) {
		cout << "Error while generating the graph" << endl;
		return ERROR;
	}


	cout << "Input:        " << input             << endl;
	cout << "Output:       " << input_coupling    << endl;
	cout << "Statistic:    " << output            << endl;
	cout << "CouplingFile: " << output_statistics << endl;
	cout << "Verbose:      " << verbose           << endl;
 	cout << "NoOuput:      " << no_output         << endl;
	 
	
	
	
	
	
	
	return SUCCESS;
}