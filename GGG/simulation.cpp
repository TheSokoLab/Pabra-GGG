//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  s i m u l a t i o n . c p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        BY TOMEK    
//        Oct 2011
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include"simulation.hpp"

using namespace std;



// ***********
// CONSTRUCTOR
// ***********
simulation::simulation(run_parameter* input){   
    
    // Copy no. of observables and dimensions
    no_obs    = (*input).no_obs;
    no_r_nucl = (*input).no_r_nucl;
    no_z_nucl = (*input).no_z_nucl;
    subv_size = (*input).subv_size;
    
    interrupts = (*input).interrupts;
    
    // Construct the double arrays of random_variable classes (1. index = species, 2. index = position)
    obs_config = new random_var*[no_obs];
    avg_obs_config = new random_var*[no_obs];
    for(unsigned int s=0; s<no_obs; s++)    // for all species
    {        
        obs_config[s]       = new random_var[no_r_nucl*no_z_nucl];
        avg_obs_config[s]   = new random_var[no_r_nucl*no_z_nucl];

        if((*input).interrupts) // if there have been temporary outputs (if this is a restart) then copy existing measurements
        {
            obs_config[s]       = (*input).obs_config[s];
            avg_obs_config[s]   = (*input).avg_obs_config[s];
        }
    }
    
    // Copy boundary thresholds
    threshold = new double[no_obs];
    for(unsigned int s=0; s<no_obs; s++)    // for all species
    {
        threshold[s] = (*input).obs[s].threshold;
    }
    
    // Construct center-of-mass measurement arrays
    obs_int = new double[no_obs];
    CoM     = new double[no_obs];
       
    // Construct boundary position histograms
    // We assume 2 boundaries per species
    obs_boundary_pos = new histo1D[2*no_obs];
    obs_boundary_no = new histo1D[2*no_obs];
    for(unsigned int s=0; s<2*no_obs; s++)  // for all species
    {
        obs_boundary_pos[s].init(1,no_z_nucl-1);
        obs_boundary_no[s].init(1,no_z_nucl-1);
        
        if((*input).interrupts){  // if there have been temporary outputs (if this is a restart) then copy existing measurements

            for(unsigned int i=1;i<no_z_nucl;i++){ // go through the whole histogram

            obs_boundary_pos[s].events_to( i, (*input).obs_boundary_pos[s].events_at(i) );
            obs_boundary_pos[s].weight_to( i, (*input).obs_boundary_pos[s].weight_at(i) );
            obs_boundary_no[s].events_to( i, (*input).obs_boundary_no[s].events_at(i) );
            obs_boundary_no[s].weight_to( i, (*input).obs_boundary_no[s].weight_at(i) );
            }
        }
    }    
    // Construct array for measured number of boundary crossings, needed in the production run do-loop later
    no_obs_boundary = new unsigned int[2*no_obs];
    
}


// ***********
// DESTRUCTOR
// ***********
simulation::~simulation(){
    
    // Clean up memory
    
    // 2D dynamic arrays
    for(unsigned int s=0; s<no_obs; s++)    // for all species
    { 
        delete[] obs_config[s];
        delete[] avg_obs_config[s];
    }
    delete[] obs_config;
    delete[] avg_obs_config;
    
    // 1D dynamic arrays
    delete[] threshold;
    
    delete[] obs_boundary_pos;
    delete[] obs_boundary_no;
    delete[] no_obs_boundary;
    
    delete[] obs_int;
    delete[] CoM;
    
    if(traf)
    {      
        traf.flush();
        traf.close();
    }
    
}


// **************
// INITIALIZATION
// **************
void simulation::init(run_parameter* input){
  
    // Construct and initialize the simulated system
    mtrand.seed((*input).random_seed);
    fly.init_cortex(mtrand, (*input));
    
    // Calculate break_steps
    break_steps = ((*input).rlx_steps + (*input).main_steps) / (*input).max_interrupts;
    // And no. of steps between instant outputs
    instant_steps = max((unsigned long)1, (unsigned long)(( (*input).main_steps) / (*input).instant_outputs ));
    
    // Now check if there have been previous temporary outputs of the
    // simulation. If this is the case, this simulation is a restart.
    // Then the random variables, the state of the system, the state
    // of the random number generator and the next-event times tree
    // have to be initialized with the values of the last output.
    
    // Initialize times and random number generator
    if((*input).interrupts)
    {              
          // Restore next event times tree
          fly.nuclear_queue.set_size((*input).tree_size);
          for(int n=0; n<(*input).tree_size; n++){

              fly.nuclear_queue.set_value(n, (*input).tree_value[n]);
              fly.nuclear_queue.set_index(n, (*input).tree_index[n]);
          }
          for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
          {
              for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
              {
                  fly.set_time(r_nucl, z_nucl, (*input).event_time[z_nucl+r_nucl*no_z_nucl]);
              }
          }
          
          // Restore random number generator state
          mtrand.load((*input).generator_state);
          
          // Restore time variables
          ac_run_step = (*input).ac_run_step;
          ac_run_time = (*input).ac_run_time;
          last_time = (*input).last_time;
          next_time = (*input).next_time;
          
          ac_meas_step = (*input).ac_meas_step;
          meas_time = (*input).meas_time;
          bgn_meas = (*input).bgn_meas;
          end_meas = (*input).end_meas;          
          
          // Check whether the system was relaxed is done later
          // to make sure the last update of the last output
          // is not measured again        
    }
    else    // Case without previous interrupts
    {
          ac_run_step = 0;
          ac_run_time = 0.0;
          last_time = 0.0;
          next_time = fly.nuclear_queue.get_value(0);
          
          ac_meas_step = 0;
          meas_time = (*input).meas_time;
          bgn_meas = 0.0;          
          end_meas = (*input).meas_time;
    }

    // Set counter for this simulations run to zero and flags to false
    // and then start the production and averaging loop
    if( (*input).rlx_steps == 0 || (*input).rlx_time == 0.0 )
    {
        // If the relaxation time is set to zero (no relaxation required)
        // an instant output is generated before the first simulation step
        // and the "relaxed" flag is set from the very beginning
        relaxed_flag  = true;
        instant_flag  = true;      
    }
    else
    {
        // Default case
        relaxed_flag  = false;    
        instant_flag  = false;
    }
    
    measure_flag = false;
    this_sim_run_step = 0;   

    // Decide whether to open the observables trajectory file
    if(mode!=AVERAGES && mode!=FFS)
    {
        // Construct stream for trajectory file
        (*input).open_output(traf, (*input).traject_file, ios_base::app, 7);
    }
  
}

// **************
// SIMULATION RUN
// **************
void simulation::run(run_parameter* input)
{
  
  // This is the main simulation loop
  //
  // At every update it is checked whether:
  //
  //    (1) a measurement should be performed for data acquisition
  //    (2) an instantaneous output should be generated
  //    (3) the simulation should be stopped
  //
  // making use of the flags defined properly by the init() routine.
    
  do {
            
      // Before any update check whether measurement or output has to be done
      if(relaxed_flag)
      {
          // *** MEASURMENT STEP ***
          if(measure_flag && mode!=TRAJECTORY && mode!=FFS)     measure();
          
          // *** INSTANTANEOUS OUTPUT ***
          // Do this independently of averaging given that the system is relaxed
          if(instant_flag && mode!=AVERAGES && mode!=FFS)       output_instant();
              
      }
      
      // *** BREAK CRITERION ***
      // Break directly after measurement to make sure the last update gets measured
      if(this_sim_run_step == break_steps)    break;
      
      // *******************************
      // *** GILLESPIE-STEP / UPDATE ***
      // *******************************
      // This part is done both in the relaxation and measurement run
      fly.gillespie_step(mtrand);

      last_time = next_time;
      ac_run_time = last_time;
      next_time = fly.nuclear_queue.get_value(0);

      this_sim_run_step++;
      ac_run_step++;
      ac_meas_step++;      
      
      // *** CHECK FOR ACTIONS TO DO BEFORE NEXT UPDATE ***
      // Check whether to start measurement from the next loop on
      if(!relaxed_flag && (ac_run_step >= (*input).rlx_steps || ac_run_time >= (*input).rlx_time) )
      {
        
        relaxed_flag  = true;
        
        // Averaging marks are set only after relaxation, no measurement desired before
        // VERY IMPORTANT: This has to be done to get the first weight_time right!
        // However only set these initial values only ONCE after relaxation!
        // Therefore check whether bgn_meas is still 0.0 as initially before relaxation:              
        if(bgn_meas == 0.0){
          
              // Initialize the averaging marks
              bgn_meas = last_time;
              end_meas = last_time + (*input).meas_time;
              ac_meas_step = 0;
        }              
      }
      
      // Check whether the next loop is a measurement loop (then: measure and restart averaging)
      if( relaxed_flag && (ac_meas_step >= (*input).meas_steps || ac_run_time >= end_meas) )
        
              measure_flag = true;
      
      // Check whether in the next loop instantaneous data should be written
      // Flag output_all discriminates between output at every averaging step or only (*input).instant_outputs times per sim.
      if( (*input).output_all &&
            (   measure_flag 
            ||  ac_run_step == (*input).rlx_steps
            )
        )          
              instant_flag = true;
      
      if( not((*input).output_all) &&
            (  (ac_run_step-(*input).rlx_steps)%instant_steps == 0         // instant output every instant_steps
            ||  ac_run_step == (*input).rlx_steps                          // ..and directly after relaxation
            )
          )
              instant_flag = true;
              // (note that instant=true only has an effect if relaxed=true)

    
    } while(next_time <= (*input).main_time + (*input).rlx_time);
    // end of main simulation loop    
}

// *********************************
// OBSERVABLE THRESHOLDS CALCULATION
// *********************************
void simulation::calculate_thresholds(run_parameter* input)
{
    // This routine typically is run at the end of the relaxation run
    // to determine the thresholds for the boundary measurement based
    // on the system-wide max. of the observables
    
    double obs_max;
    double obs_avg;    
    
    for(unsigned int s=0; s<no_obs; s++){
        
        obs_max = 0.0;  // Initialize max for this observable
      
        // Calculate the max. of the angular average
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
            
            // First calulate the angular average for this z-value
            obs_avg = 0.0;            
            for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
        
                obs_avg += (double)fly.get_observable(r_nucl,z_nucl,s);          
            }
            obs_avg = obs_avg / no_r_nucl;
            
            // Update max.
            if(obs_avg > obs_max)    obs_max = obs_avg;
              
        }   // z
    
        // Set the threshold for this observable to be 1/2 of the max.
        threshold[s] = 0.5*obs_max;
        // Also change the values in the master input parameter class
        // to ensure that the right infoline is written at the end
        (*input).threshold[s] = 0.5*obs_max;
        
    }   // s           
    
}

// **********************
// OBSERVABLE MEASUREMENT
// **********************
void simulation::measure(void)
{
    // Calculate the statistical weight of this measurement
    weight_time = last_time - bgn_meas;

    // (1.) Add configuration to random variable
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
    for(unsigned int s=0; s<no_obs; s++)  // for all species
    {
        obs_config[s][z_nucl+r_nucl*no_z_nucl].rnd_add(double(fly.get_observable(r_nucl,z_nucl,s)), weight_time);
        avg_obs_config[s][z_nucl].rnd_add(double(fly.get_observable(r_nucl,z_nucl,s)), weight_time);            
    }
    }
    }

    // (2.) Find position and number of threshold crossings
    // First find concentration "center of mass"
    for(unsigned int s=0; s<no_obs; s++)
    {
        obs_int[s] = 0;
        CoM[s]   = 0;
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
        for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
            
            obs_int[s] += (double)fly.get_observable(r_nucl,z_nucl,s);
            CoM[s]     += (double)(fly.get_observable(r_nucl,z_nucl,s) * z_nucl);
        }
        }
        
        CoM[s] = CoM[s]/obs_int[s];
    }    
    // Reset measurement array
    for(unsigned int s=0; s<2*no_obs; s++) no_obs_boundary[s] = 0;    
    // Measure it   
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        double is_this_above;
        double is_next_above;
        double is_boundary;
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl-1; z_nucl++){
        for(unsigned int s=0; s<no_obs; s++)  // for all species
        {
            // Check whether the current nucleus and/or the next one are above the specified threshold
            is_this_above = fly.get_observable(r_nucl,z_nucl,s) - threshold[s];
            is_next_above = fly.get_observable(r_nucl,z_nucl+1,s) - threshold[s];
            
            // If one of them is above and the other is not, we have a boundary here
            if(is_this_above*is_next_above < 0)
            {
                // A boundary has been found
                // Now determine whether it's the anterior or posterior
                if( 1.0*z_nucl+0.5 <= 1.0*CoM[s])
                {
                    // anterior boundary (before CoM)
                    obs_boundary_pos[s].add(z_nucl+1, weight_time);
                    no_obs_boundary[s]++;
                }
                else
                {
                    // posterior boundary (beyond CoM)
                    obs_boundary_pos[no_obs+s].add(z_nucl+1, weight_time);
                    no_obs_boundary[no_obs+s]++;
                }
              }
        }
        }
    }
    // Add to histograms
    for(unsigned int s=0; s<2*no_obs; s++)    obs_boundary_no[s].add(no_obs_boundary[s], weight_time);

    // (3.) Reset measurement initial condition
    measure_flag = false;
    ac_meas_step = 0;
    bgn_meas = last_time;
    end_meas = last_time + meas_time;

    // Count total number of measurements in this simulation
    N_measurements++;
    
}

// *******************************
// REACTION COORDINATE CALCULATION
// *******************************
double* simulation::calculate_reaction_coordinate(void)
{
    // This is still hard-coded to avoid (for now) implementation
    // of a function parser; might be generalized at some point
    
    enum Species { Hb_tot   = 2,
                   Kr_tot   = 7,
                   Kni_tot  = 12,
                   Gt_tot   = 17,
                   Hb_prom  = 3, 
                   Kr_prom  = 8,
                   Kni_prom = 13,
                   Gt_prom  = 18  };      
    
    double  int_tot_sum;
    double  int_prom_sum;
    int     rc_no = 12 + no_z_nucl; // number of reaction coordinates
    double* rc_array = new double[rc_no+1]; // reaction coordinates return array   
    
    double* av_temp  = new double[no_z_nucl+1]; // r-averages temporary array
    double* av_store = new double[no_z_nucl+1]; // r-averages storage
    double  base;
    
    // Initialize cumulative copy number array
    for(unsigned int s=0; s<no_obs; s++)    obs_int[s] = 0;
    
    // Integrate over the whole system
    for(unsigned int s=0; s<no_obs; s++){
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
    
        obs_int[s] += (double)fly.get_observable(r_nucl,z_nucl,s);
      
    }   // z
    }   // r
    }   // s
    
    // Calculate the total sums over the whole system
    int_tot_sum  = (double)(obs_int[Hb_tot] + obs_int[Kr_tot] + obs_int[Kni_tot] + obs_int[Gt_tot]);
    int_prom_sum = (double)(obs_int[Hb_prom] + obs_int[Kr_prom] + obs_int[Kni_prom] + obs_int[Gt_prom]);
    
    // Calculate and store the averages of conc. around the circumference
    for(unsigned int s=0; s<no_obs; s++){
      
      // First reset averaging array
      for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        
          av_temp[z_nucl] = 0.0;
        
      // Now for each z-pos. ...
      for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
        
          // ... do the averaging
          for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
          
            av_temp[z_nucl] += (double)fly.get_observable(r_nucl,z_nucl,s);
            
          }   // r          
          
          // Construct a 16-digit number that contains all 
          // relevant average copy numbers at this z
          switch(s){
          
            case Hb_tot :
              
              base = 1.0e0;
              break;
              
            case Kr_tot :
              
              base = 1.0e3;
              break;
            
            case Kni_tot :
              
              base = 1.0e6;
              break;
              
            case Gt_tot :
              
              base = 1.0e9;
              break;
              
            default :
              
              base = 0.0;
              break;
              
          } // switch s
              
          // Add to the constructed "word" for this z-position
          av_store[z_nucl] += (double)( base * (int)(1.0/no_r_nucl*av_temp[z_nucl]) );
          
      } // z      
    }   // s
    
    // Construct the RC array
    rc_array[0] = rc_no;
    rc_array[1] = obs_int[Hb_tot];
    rc_array[2] = obs_int[Kr_tot];
    rc_array[3] = obs_int[Kni_tot];
    rc_array[4] = obs_int[Gt_tot];
    rc_array[5] = int_tot_sum;
    rc_array[6] = (obs_int[Hb_tot] + max(obs_int[Kr_tot], obs_int[Gt_tot])) / int_tot_sum;
    
    rc_array[7]  = obs_int[Hb_prom];
    rc_array[8]  = obs_int[Kr_prom];
    rc_array[9]  = obs_int[Kni_prom];
    rc_array[10] = obs_int[Gt_prom];
    rc_array[11] = int_prom_sum;
    rc_array[12] = (obs_int[Hb_prom] + max(obs_int[Kr_prom], obs_int[Gt_prom])) / int_prom_sum;

    // Copy averages array to RC array
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++){
      
      rc_array[13+z_nucl] = av_store[z_nucl];
      
    }
      
    return rc_array;
    
}

// *************************************
// OUTPUT OF INSTANTANEOUS CONFIGURATION
// *************************************
void simulation::output_instant(void)
{
  
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        // Calculate relative r-position and r-distance
        double r_pos = double(r_nucl+0.5)/double(no_r_nucl);
        double r_dst = double(r_nucl+0.5)*double(subv_size);
        
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            // Calculate relative z-position and z-distance
            double z_pos = double(z_nucl+0.5)/double(no_z_nucl);
            double z_dst = double(z_nucl+0.5)*double(subv_size);
            
            // Output position and time
            traf << subv_size << " "
                 << no_r_nucl << " " << r_nucl << " " << r_pos << " " << r_dst << " "
                 << no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " "
                 << ac_run_step << " " << interrupts*break_steps + this_sim_run_step << " "
                 << ac_run_time << " " << 0.5*(next_time+last_time) << " " << last_time << " ";
                  
            // Output species copy numbers (starting at column 14, counting from 0)
            for(unsigned int s=0; s<no_obs; s++)   // for all species
            {
                traf << double(fly.get_observable(r_nucl,z_nucl,s)) << " ";
            }
            traf << endl;
            
        }   // end z-loop
    }       // end r-loop
              
    // Reset instant flag to false to avoid output in the next loop
    instant_flag = false;
    
}

// **********************************
// OUTPUT OF COMPLETE TEMPORARY STATE
// **********************************
void simulation::output_temp_file(run_parameter* input)
{
    // The file created by this routine contains the complete
    // information to restart the simulation at precisely the
    // state at which it was interrupted.
    
    // Store random number generator state
    MTRand::uint32 generator_state[mtrand.SAVE];
    mtrand.save(generator_state);
    
    // Open temporary output file
    (*input).open_output(tempf, (*input).temp_file, 0, 10);
    
    // Increment interrupts number and start writing to temporary file
    (*input).interrupts++;
    tempf << (*input).interrupts << endl;
    tempf << (*input).random_seed << endl;

    // Save mersenne twister random number generator state
    tempf << mtrand.SAVE << endl;
    for(unsigned int i=0; i<mtrand.SAVE; i++)   tempf << generator_state[i] << " ";
    tempf << endl;

    for(unsigned int s=0; s<no_obs; s++)   // for all (measured) species
    {
        // Save random variables
        for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {              
            // Non-averaged observables
            tempf << obs_config[s][z_nucl+r_nucl*no_z_nucl].events << " "
                  << obs_config[s][z_nucl+r_nucl*no_z_nucl].weight << " "
                  << obs_config[s][z_nucl+r_nucl*no_z_nucl].max_moments;
            // also save moments as those are calculated on the fly
            for(int j = 0; j<obs_config[s][z_nucl+r_nucl*no_z_nucl].max_moments; j++)
                  tempf << " " << obs_config[s][z_nucl+r_nucl*no_z_nucl].moments[j];
            tempf << endl;

            // r-averaged observables
            tempf << avg_obs_config[s][z_nucl].events << " "
                  << avg_obs_config[s][z_nucl].weight << " " 
                  << avg_obs_config[s][z_nucl].max_moments;
            // also save moments as those are calculated on the fly
            for(int j = 0; j<avg_obs_config[s][z_nucl].max_moments; j++)
                  tempf << " " << avg_obs_config[s][z_nucl].moments[j];
            tempf << endl;
        }
        }

        // Save boundary histogram counts and weights
        for(unsigned int i=1; i<no_z_nucl; i++)
        {
            tempf << obs_boundary_pos[s].events_at(i) << " " << obs_boundary_pos[s].weight_at(i) << " "
                  << obs_boundary_no[s].events_at(i)  << " " << obs_boundary_no[s].weight_at(i)  << " "
                  << obs_boundary_pos[no_obs+s].events_at(i) << " " << obs_boundary_pos[no_obs+s].weight_at(i) << " "
                  << obs_boundary_no[no_obs+s].events_at(i)  << " " << obs_boundary_no[no_obs+s].weight_at(i)
                  << endl;
        }
    } // species loop (s)

    // Save current nuclei configurations and next-event times
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++){
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
    {
        // for all (measured and unmeasured) species
        for(unsigned int spec = 0; spec < (*input).no_spec; spec++)
            tempf << fly.get_config(r_nucl,z_nucl,spec) << " ";

        tempf << fly.get_time(r_nucl,z_nucl) << endl;
    }            
        tempf << endl;
    }
    
    // Save binary tree of next event times
    tempf << ac_run_step << " " << ac_run_time << " " << last_time << " " << next_time << endl;
    if((*input).meas_time > 0.0 || (*input).meas_steps > 0) // write this line only if averaging was on
          tempf << ac_meas_step << " " << bgn_meas << " " << end_meas << endl;
    tempf << fly.nuclear_queue.get_size() << endl;
    for(unsigned int n = 0; n < fly.nuclear_queue.get_size(); n++)
          tempf << fly.nuclear_queue.get_value(n) << " " << fly.nuclear_queue.get_index(n) << endl;
    
    tempf << N_measurements << endl;

    // Finalize the filestream
    tempf.flush();
    (*input).close_output(tempf);
    
    // Set name of "current_temp_file"; at the end the temp_file is backuped (copied) to this name
    strcpy((*input).current_temp_file, (*input).temp_file);
    strcat((*input).current_temp_file, ".");
    sprintf(tmp_string, "%012lu", ac_run_step);
    strcat((*input).current_temp_file, tmp_string);
    
    // Construct copy command (UNIX)
    strcpy(tmp_string, "cp ");
    strcat(tmp_string, (*input).temp_file);
    strcat(tmp_string, " ");
    strcat(tmp_string, (*input).current_temp_file);
    // Execute it (this will work only on UNIX systems)
    int syserr;        
    syserr=system(tmp_string);
    
}

// *********************************************
// OUTPUT OF INPUT FILE TO BE USED FOR A RESTART
// *********************************************
void simulation::output_input_file(run_parameter* input)
{ 
    // Open the file
    ofstream new_input_file((*input).new_input_file, ios::out);
    new_input_file.precision(5);

    // Output general parameters
    new_input_file << (*input).what_to_do       << endl;
    new_input_file << (*input).max_interrupts   << endl;
    new_input_file << (*input).init_random_seed << endl;

    new_input_file << (*input).rlx_time    << endl;
    new_input_file << (*input).main_time   << endl;
    new_input_file << (*input).meas_time   << endl;
    new_input_file << (*input).rlx_steps   << endl;
    new_input_file << (*input).main_steps  << endl;
    new_input_file << (*input).meas_steps  << endl;
    new_input_file << (*input).no_trajec   << endl;

    new_input_file << (*input).subv_size << endl;
    new_input_file << (*input).no_r_nucl << endl;
    new_input_file << (*input).no_z_nucl << endl;

    new_input_file << (*input).no_reac   << endl;
    new_input_file << (*input).no_para   << endl;
    new_input_file << (*input).no_spec   << endl;
    new_input_file << (*input).no_obs    << endl;
    new_input_file << (*input).no_neig   << endl;
    
    // Thresholds
    for(unsigned int spec = 0; spec < (*input).no_spec; spec++)
    {
        new_input_file << (*input).threshold[spec] << endl;
    }
    
    // Reaction rates
    for(unsigned int r_nucl = 0; r_nucl < (*input).no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < (*input).no_z_nucl; z_nucl++)
        {
            for(unsigned int para = 0; para < (*input).no_para; para++)
            {
                new_input_file << (*input).reac_parameter[para+z_nucl*(*input).no_para+r_nucl*(*input).no_z_nucl*(*input).no_para] << " ";
            }
            
            new_input_file << endl;
        }
    }
    
    // Diffusion rates
    for(unsigned int r_nucl = 0; r_nucl < (*input).no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < (*input).no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < (*input).no_spec; spec++)
            {
                new_input_file << (*input).diff_rate_const[spec+z_nucl*(*input).no_spec+r_nucl*(*input).no_z_nucl*(*input).no_spec] << " ";
            }
            
            new_input_file << endl;
        }
    }
    
    // Configuration
    for(unsigned int r_nucl = 0; r_nucl < (*input).no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < (*input).no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < (*input).no_spec; spec++)
            {
                // This is the part that differs from the intial input file
                new_input_file << fly.get_config(r_nucl,z_nucl,spec) << " ";
            }
            
            new_input_file << endl;
        }
    }
    
    new_input_file << (*input).output_all << endl;
    new_input_file << (*input).instant_outputs << endl;
    
    // Clean up
    new_input_file.flush();
    new_input_file.close();
    
}

// ***********************************
// AVERAGING OF OBSERVABLES AND OUTPUT
// ***********************************
void simulation::output_averages(run_parameter* input)
{
    // Open output files
    ofstream average;       (*input).open_output(average, (*input).average_file, 1);
    ofstream bound_dist;    (*input).open_output(bound_dist, (*input).bnd_dst_file, 1);

    // Define measurement arrays
    double** ttl_obs_config = new double*[no_obs];
    double** sum_obs_config = new double*[no_obs];
    double** prc_obs_config = new double*[no_obs];
    for(unsigned int s=0; s<no_obs; s++)    // for all species
    {
        ttl_obs_config[s] = new double[no_r_nucl];
        sum_obs_config[s] = new double[no_r_nucl];
        prc_obs_config[s] = new double[no_r_nucl];
    }
    
    
    // *** SPATIAL PROFILES ***
    // Write info line
    average << "### INFO LINE: DATA FORMAT #########################################" << endl;
    average << "### First 9 rows are: " << endl;
    average << "### subv_size | " << "no_r_nucl | " << "r_nucl | " << "r_pos | " << "r_dst | ";
    average <<     "no_z_nucl | " << "z_nucl | " << "z_pos | " << "z_dst | ";
    average << endl;
    average << "### Next rows are the following 3 quantities for all observables: " << endl;
    average << "### obs_config[s][z_nucl+r_nucl*no_z_nucl].mean() | "
            << "obs_config[s][z_nucl+r_nucl*no_z_nucl].variance() | "
            << "prc_obs_config[s][r_nucl] | " << endl;
    average << "### Observable list: " << endl;
    for(int s=0; s<no_obs; s++)
        average << "#   Observable " << s << " = " << (*input).obs[s].name << endl;
    average << "####################################################################" << endl;
    
    // *** SPATIAL PROFILES FOR EVERY "RADIUS" [r_nucl] ***
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        // First set all arrays to zero
        for(unsigned int s=0; s<no_obs; s++)    // for all species
        {
            ttl_obs_config[s][r_nucl] = 0.0;
            sum_obs_config[s][r_nucl] = 0.0;
            prc_obs_config[s][r_nucl] = 0.0;
        }
    }
    // Calculate cumulative sum for every r-position
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int s=0; s<no_obs; s++)    // for all species
                ttl_obs_config[s][r_nucl] += obs_config[s][z_nucl+r_nucl*no_z_nucl].mean();
        }
    }
    // Normalize and write fo file
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            // Calculate relative positions and distances from position indices
            double r_pos = double(r_nucl+0.5)/double(no_r_nucl);
            double z_pos = double(z_nucl+0.5)/double(no_z_nucl);
            double r_dst = double(r_nucl+0.5)*subv_size;
            double z_dst = double(z_nucl+0.5)*subv_size;

            average << subv_size << " ";
            average << no_r_nucl << " " << r_nucl << " " << r_pos << " " << r_dst << " ";
            average << no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " ";

            for(unsigned int s=0; s<no_obs; s++)    // for all species
            {
                sum_obs_config[s][r_nucl] += obs_config[s][z_nucl+r_nucl*no_z_nucl].mean();
                prc_obs_config[s][r_nucl] = sum_obs_config[s][r_nucl]/ttl_obs_config[s][r_nucl];

                average << obs_config[s][z_nucl+r_nucl*no_z_nucl].mean() << " " << obs_config[s][z_nucl+r_nucl*no_z_nucl].variance() << " " << prc_obs_config[s][r_nucl] << " ";
            }           
            average << endl;
        }
    }
    
    // Clean up
    for(unsigned int s=0; s<no_obs; s++)
    {
          delete[] ttl_obs_config[s];
          delete[] sum_obs_config[s];
          delete[] prc_obs_config[s];
    }
    delete[] ttl_obs_config;
    delete[] sum_obs_config;
    delete[] prc_obs_config;


    // *** RADIAL AVERAGES OF PROFILES ***
    // Define new measurement arrays
    double* ttl_avg_obs_config = new double[no_obs];
    double* sum_avg_obs_config = new double[no_obs];
    double* prc_avg_obs_config = new double[no_obs];

    // Write another info line
    average << "# * INFO LINE: AVERAGE DATA FORMAT **********************************************************************************************************"
            << endl;
    average << "### First 5 rows are: " << endl;
    average << "# * subv_size | " << "no_z_nucl | " << "z_nucl | " << "z_pos | " << "z_dst | " << endl;
    average << "# * Next rows are the following 3 quantities for all observables: " << endl;
    average << "# * avg_obs_config[s][z_nucl].mean() | "
            <<      "avg_obs_config[s][z_nucl].variance() | "
            <<      "avg_prc_obs_config[s] | " << endl;
    average << "# * Observable list: " << endl;
    for(int s=0; s<no_obs; s++)
        average << "# *   Observable " << s << " = " << (*input).obs[s].name << endl;
    average << "# *******************************************************************************************************************************************"
            << endl;

    // Calculate cumulative sum for every species
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
    {
        for(unsigned int s=0; s<no_obs; s++)    // for all species
            ttl_avg_obs_config[s] += avg_obs_config[s][z_nucl].mean();
    }
    // Normalize and write fo file
    for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
    {
        // Calculate relative positions and distances from position indices
        double z_pos = double(z_nucl+0.5)/double(no_z_nucl);
        double z_dst = double(z_nucl+0.5)*subv_size;

        average << "* " << subv_size << " " << no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " ";

        for(unsigned int s=0; s<no_obs; s++)    // for all species
        {
            sum_avg_obs_config[s] += avg_obs_config[s][z_nucl].mean();
            prc_avg_obs_config[s] =  sum_avg_obs_config[s]/ttl_avg_obs_config[s];
            
            average << avg_obs_config[s][z_nucl].mean() << " " << avg_obs_config[s][z_nucl].variance() << " " << prc_avg_obs_config[s] << " ";
        }
        average << endl;
    }
    
    // Clean up
    delete[] ttl_avg_obs_config;
    delete[] sum_avg_obs_config;
    delete[] prc_avg_obs_config;
    
    //Finalize
    average.flush();
    (*input).close_output(average);


    // *** BOUNDARY DISTRIBUTIONS ***
    // Info line
    bound_dist << "### INFO LINE: BOUNDARY DISTRIBUTION DATA FORMAT ################" << endl;
    bound_dist << "### First 5 rows are: " << endl;
    bound_dist << "### subv_size | " << "no_z_nucl | " << "z_nucl | " << "z_pos | "  << "z_dst | " << endl;
    bound_dist << "### Next rows are the following 3 quantities for all boundaries (first anterior, then posterior for each species): " << endl;
    bound_dist << "### obs_boundary_pos[s].average() | " << "obs_boundary_pos[s].variance() | " << "obs_boundary_pos[s].densty_at(z_nucl) | " << endl;
    bound_dist << "### Observable list: " << endl;
    for(int s=0; s<no_obs; s++)
        bound_dist << "#   Observable " << s << " = " << (*input).obs[s].name << endl;
    bound_dist << "#################################################################" << endl;    
    bound_dist << endl;

    for(unsigned int s=0; s<2*no_obs; s++)  // for all species
    {
        obs_boundary_pos[s].normalise();
        obs_boundary_no[s].normalise();
    }

    for(unsigned int z_nucl = 1; z_nucl < no_z_nucl; z_nucl++)
    {
        // Calculate relative positions and distances from position indices
        double z_pos = double(z_nucl)/double(no_z_nucl);
        double z_dst = double(z_nucl)*subv_size;

        bound_dist << subv_size << " " << no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " ";

        for(unsigned int s=0; s<no_obs; s++)    // for all species
        {
            bound_dist  << obs_boundary_pos[s].average() << " " << obs_boundary_pos[s].variance() << " " << obs_boundary_pos[s].densty_at(z_nucl) << " "
                        << obs_boundary_pos[no_obs+s].average() << " " << obs_boundary_pos[no_obs+s].variance() << " " << obs_boundary_pos[no_obs+s].densty_at(z_nucl) << " ";
                        //<< obs_boundary_no[s].average() << " " << obs_boundary_no[s].variance() << " " << obs_boundary_no[s].densty_at(z_nucl) << " ";
        }
        bound_dist << endl;
    }
    
    // Finalize
    bound_dist.flush();
    (*input).close_output(bound_dist);
    
}

