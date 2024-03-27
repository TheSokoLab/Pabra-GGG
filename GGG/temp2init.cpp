
#include<iostream>
#include "run_parameter.hpp"


using namespace std;

int main(int argc, char *argv[])
{
 
    if(argc != 5) 
    {
        cerr << "You have to pass 2 arguments to function main!" << endl;
        exit(EXIT_FAILURE);
    }
    
    
    run_parameter input;
    strcpy(input.temp_file, argv[1]);
    strcpy(input.input_file, argv[2]);
    strcpy(input.reactions_file, argv[3]);
    strcpy(input.observables_file, argv[4]);
    
    input.read_parameter();
    input.read_temp_file();
    
    cout << "# Initial conditions filtered from " << input.temp_file << " with parameters taken from " << input.input_file << endl;
    cout << "# no_r_nucl = " << input.no_r_nucl << ", no_z_nucl = " << input.no_z_nucl << ", no_spec = " << input.no_spec << endl;
    cout << "# Format: r_nucl, z_nucl, configuration[0], configuration[1], [..]" << endl;
    cout << "# " << endl;
        
    // Output configurations    
    int i=0;
    for(unsigned int r_nucl = 0; r_nucl < input.no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < input.no_z_nucl; z_nucl++)
        {
            cout << r_nucl << " " << z_nucl << " ";
            
            for(unsigned int spec = 0; spec < input.no_spec; spec++)  // for all species
                cout << input.configuration[spec+z_nucl*input.no_spec+r_nucl*input.no_z_nucl*input.no_spec] << " ";        
            cout << endl;
            
            i++;
        }
    }
    cout << "# Written " << i << " configurations." << endl;
    
     
}

