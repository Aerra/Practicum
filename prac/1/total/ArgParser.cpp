#include "ArgParser.h"

void ArgParser::PrintHelp() {
    std::cout << "Usage: "
    << "\n-h Print this information and return\n"
    << "-v if exist then print result to stdout\n"
    << "--tol=[param] relative precision of decision\n"
    << "--T=[param] Count of threads\n"
	<< "--Nx=[param] Count of Nx cell (column)\n"
    << "--Ny=[param] Count of Ny cell (raw)\n"
    << "--K1=[param] Count of quadrangle\n"
    << "--K2=[param] Count of diagonally divided quadrilateral\n"
    << "--output=[filepath] Path to file with result (if not "
	<< "found then verbose param not turn on automatically)\n"
    << "--input=[filepath] Path to file, that contain Nx, Ny, K1 "
	<< "and K2 params (format [Nx]\\n[Ny]\\n[K1]\\n[K2]\\n[T]\\n[tol])\n"
	<< "-d if exist enable debug mode\n";
}

// Use ONLY for incoming arguments parsing
std::string getCmdOption(int argc, char* argv[], const std::string& option, \
						bool &find)
{   
    std::string cmd;
    find = false;
    for( int i = 0; i < argc; ++i)
    {
        std::string arg = argv[i];
        if(0 == arg.find(option))
        {
            find = true;                
            cmd = arg.substr(option.length());
            return cmd;
        }
    }
    return cmd;
}

// Function that parses incoming arguments and saves them in class ArgParser
// return false if main function should exit after calling Parse()
// return true in other cases
bool ArgParser::Parse(int argc, char *argv[]) {
	bool find = false;
	std::string help_s = getCmdOption(argc, argv, "-h", find);
	if (find) {
		PrintHelp();
	    return false;
	}
	
	std::string verbose_s = getCmdOption(argc, argv, "-v", find);
	if (find) {
	    verbose = true;
	}

	std::string debug_s = getCmdOption(argc, argv, "-d", find);
	if (find) {
	    debug = true;
	}
	
	input_file = getCmdOption(argc, argv, "--input=", find);	
	// no input file - so read from argv
	if (!find) {
		std::string Nx_s = getCmdOption(argc, argv, "--Nx=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter Nx!\n";
		    return false;
		}
		Nx = atoi(Nx_s.c_str());
		if (Nx == 0) {
		    std::cout << "Invalid Nx has been founded!\n";
		    return false;
		}				

		std::string Ny_s = getCmdOption(argc, argv, "--Ny=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter Ny!\n";
		    return false;
		}
		Ny = atoi(Ny_s.c_str());
		if (Ny == 0) {
		    std::cout << "Invalid Ny has been founded!\n";
		    return false;
		}				

		std::string K1_s = getCmdOption(argc, argv, "--K1=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter K1!\n";
		    return false;
		}
		K1 = atoi(K1_s.c_str());
		if (K1 == 0) {
		    std::cout << "Invalid K1 has been founded!\n";
		    return false;
		}				

		std::string K2_s = getCmdOption(argc, argv, "--K2=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter K2!\n";
		    return false;
		}
		K2 = atoi(K2_s.c_str());
		if (K2 == 0) {
		    std::cout << "Invalid K2 has been founded!\n";
		    return false;
		}				

		std::string T_s = getCmdOption(argc, argv, "--T=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter T!\n";
		    return false;
		}
		T = atoi(T_s.c_str());
		if (T == 0) {
		    std::cout << "Invalid T has been founded!\n";
		    return false;
		}				

		std::string tol_s = getCmdOption(argc, argv, "--tol=", find);
		if (!find) {
			std::cout << "Can't find required parameter required parameter tol!\n";
		    return false;
		}
		tol = atof(tol_s.c_str());
		if (tol <= 0) {
		    std::cout << "Invalid tol has been founded!\n";
		    return false;
		}				
	} else {
		bool ok = true;

		// Read input file
		setlocale(LC_ALL, "rus");
		std::ifstream in_file (input_file.c_str());
		if (!in_file) {
		    std::cout << "Can't open input file " << input_file << "\n";
			return false;
		}

		std::string Nx_s;
		if (getline(in_file, Nx_s)) {
			Nx = atoi(Nx_s.c_str());
			if (Nx == 0) {
			    std::cout << "Invalid Nx has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter Nx!\n";
			ok = false;
		}

		std::string Ny_s;
		if (getline(in_file, Ny_s)) {
			Ny = atoi(Ny_s.c_str());
			if (Ny == 0) {
			    std::cout << "Invalid Ny has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter Ny!\n";
			ok = false;
		}

		std::string K1_s;
		if (getline(in_file, K1_s)) {
			K1 = atoi(K1_s.c_str());
			if (K1 == 0) {
			    std::cout << "Invalid K1 has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter K1!\n";
			ok = false;
		}

		std::string K2_s;
		if(getline(in_file, K2_s)) {
			K2 = atoi(K2_s.c_str());
			if (K2 == 0) {
			    std::cout << "Invalid K2 has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter K2!\n";
			ok = false;
		}

		std::string T_s;
		if (getline(in_file, T_s)) {
			T = atoi(T_s.c_str());
			if (T == 0) {
			    std::cout << "Invalid T has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter T!\n";
			ok = false;
		}

		std::string tol_s;
		if (getline(in_file, tol_s)) {
			tol = atof(tol_s.c_str());
			if (tol <= 0) {
			    std::cout << "Invalid tol has been founded!\n";
				ok = false;
			}				
		} else {
		    std::cout << "Can't find required parameter tol!\n";
			ok = false;
		}

		in_file.close();
		if (!ok) { return false; }
	}

	output_file = getCmdOption(argc, argv, "--output=", find);
	if (find) {
	    setlocale(LC_ALL, "rus");
	    std::ofstream out_file (output_file.c_str());
	    if (!out_file) {
	        std::cout << "Can't open output file " << output_file << "\n";
	        return false;
	    }
	    out_file.close();
	}

	return true;
}
