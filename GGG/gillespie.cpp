//  ~~~~~~~~~~~~~~~~~~~~~~~~~
//  g i l l e s p i e . c p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~

#include <iostream>
//#include "limits.h"

#include "parameters.hpp"
#include "simulation.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2) 
    {
        cerr << "You have to pass arguments to function main!" << endl;
        exit(EXIT_FAILURE);
    }
    
    run_parameter input;    
    input.read_filenames(argc, argv);
    input.read_parameter();
    input.read_temp_file();
    
    simulation sim(&input);

    switch(input.what_to_do)
    {
        case 0:
            //trajectory(input);  // OLD VERSION    
            sim.mode = TRAJECTORY;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            break;
            
        case 1:
            //boundary(input);    // OLD VERSION
            sim.mode = AVERAGES;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            // If this is the last output, run the averaging
            // Note that input.max_interrupts initially contains the number of previous outputs
            // This number is increased by sim.output_temp_file
            if(input.interrupts>=input.max_interrupts){
              
                sim.calculate_thresholds(&input);  
                sim.output_averages(&input);
            }                
            
            break;
            
        case 2:
            //profile(input);     // OLD VERSION           
            sim.mode = FFS;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            

            cout << sim.last_time                       << " ";
            cout.precision(16);
            // Output the reaction coordinates; first entry of array stores their number
            for(int i=1; i<=sim.calculate_reaction_coordinate()[0]; i++)
            {
                cout << sim.calculate_reaction_coordinate()[i] << " ";
            }
            cout << endl;
            cout.flush();
            
            break;
            
        case 3:
            //allinfo(input);   // OLD VERSION
            sim.mode = ALL;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            // If this is the last output, run the averaging
            // Note that input.max_interrupts initially contains the number of previous outputs
            // This number is increased by sim.output_temp_file
            if(input.interrupts>=input.max_interrupts){
              
                sim.calculate_thresholds(&input);  
                sim.output_averages(&input);
            }
            
            cout << "lambda = " << sim.calculate_reaction_coordinate()[4] << endl;
            cout.flush();
            
            break;

        default :
            cerr << "Don't know what_to_do. Select [0-3]!\n";
            
    }
    return(EXIT_SUCCESS);
}
