#include "ArgParser.h"

void ArgParser::PrintHelp() {
    std::cout << "Usage: " << "\n";
    std::cout << "-h Print this information and return\n";
    std::cout << "-v if exist then print result to stdout\n";
    std::cout << "--tol=[param] relative precision of decision\n";
    std::cout << "--T=[param] Count of threads\n"; std::cout << "--Nx=[param] Count of Nx cell (column)\n";
    std::cout << "--Ny=[param] Count of Ny cell (raw)\n";
    std::cout << "--K1=[param] Count of quadrangle\n";
    std::cout << "--K2=[param] Count of diagonally divided quadrilateral\n";
    std::cout << "--output=[filepath] Path to file with result (if not " <<
			"found then verbose param not turn on automatically)\n";
    std::cout << "--input=[filepath] Path to file, that contain Nx, Ny, K1 " <<
			"and K2 params (format [Nx]\\n[Ny]\\n[K1]\\n[K2]\\n[T]\\n[tol])\n";
}

// ONLY for parse incoming arguments
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
            //std::size_t found = arg.find_last_of(option);
            cmd = arg.substr(option.length());
            return cmd;
        }
    }
    return cmd;
}

// Function that parse incoming arguments and save them in class ArgParser
// return false if after Parse() calling function should exit
// return true in other case
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
	
	input_file = getCmdOption(argc, argv, "--input=", find);	
	if (!find) {
		std::string Nx_s = getCmdOption(argc, argv, "--Nx=", find);
		if (!find) { std::cout << "Can't find Nx , required!\n";
		    return false;
		}
		Nx = atoi(Nx_s.c_str());
		if (Nx == 0) {
		    std::cout << "Find invalid Nx!\n";
		    return false;
		}				

		std::string Ny_s = getCmdOption(argc, argv, "--Ny=", find);
		if (!find) {
		    std::cout << "Can't find Ny, required!\n";
		    return false;
		}
		Ny = atoi(Ny_s.c_str());
		if (Ny == 0) {
		    std::cout << "Find invalid Ny!\n";
		    return false;
		}				

		std::string K1_s = getCmdOption(argc, argv, "--K1=", find);
		if (!find) {
		    std::cout << "Can't find K1, required!\n";
		    return false;
		}
		K1 = atoi(K1_s.c_str());
		if (K1 == 0) {
		    std::cout << "Find invalid K1!\n";
		    return false;
		}				

		std::string K2_s = getCmdOption(argc, argv, "--K2=", find);
		if (!find) {
		    std::cout << "Can't find K2, required!\n";
		    return false;
		}
		K2 = atoi(K2_s.c_str());
		if (K2 == 0) {
		    std::cout << "Find invalid K2!\n";
		    return false;
		}				

		std::string T_s = getCmdOption(argc, argv, "--T=", find);
		if (!find) {
		    std::cout << "Can't find T, required!\n";
		    return false;
		}
		T = atoi(T_s.c_str());
		if (T == 0) {
		    std::cout << "Find invalid T!\n";
		    return false;
		}				

		std::string tol_s = getCmdOption(argc, argv, "--tol=", find);
		if (!find) {
		    std::cout << "Can't find tol, required!\n";
		    return false;
		}
		tol = atof(tol_s.c_str());
		if (tol == 0) {
		    std::cout << "Find invalid tol!\n";
		    return false;
		}				
	} else {
		bool not_ok = false;

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
			    std::cout << "Find invalid Nx!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find Nx , required!\n";
			not_ok = true;
		}

		std::string Ny_s;
		if (getline(in_file, Ny_s)) {
			Ny = atoi(Ny_s.c_str());
			if (Ny == 0) {
			    std::cout << "Find invalid Ny!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find Ny, required!\n";
			not_ok = true;
		}

		std::string K1_s;
		if (getline(in_file, K1_s)) {
			K1 = atoi(K1_s.c_str());
			if (K1 == 0) {
			    std::cout << "Find invalid K1!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find K1, required!\n";
			not_ok = true;
		}

		std::string K2_s;
		if(getline(in_file, K2_s)) {
			K2 = atoi(K2_s.c_str());
			if (K2 == 0) {
			    std::cout << "Find invalid K2!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find K2, required!\n";
			not_ok = true;
		}

		std::string T_s;
		if (getline(in_file, T_s)) {
			T = atoi(T_s.c_str());
			if (T == 0) {
			    std::cout << "Find invalid T!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find T, required!\n";
			not_ok = true;
		}

		std::string tol_s;
		if (getline(in_file, tol_s)) {
			tol = atof(tol_s.c_str());
			if (tol == 0) {
			    std::cout << "Find invalid tol!\n";
				not_ok = true;
			}				
		} else {
		    std::cout << "Can't find tol, required!\n";
			not_ok = true;
		}

		in_file.close();
		if (not_ok) { return false; }	
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
