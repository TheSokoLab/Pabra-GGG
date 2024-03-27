#ifndef _RUN_PARAMETER_HPP
#define _RUN_PARAMETER_HPP

#include<ctime>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<assert.h>

#include"./Tools/random_variable.hpp"
#include"./Tools/histo1D.hpp"
#include"./Tools/mersenne_twister.hpp"
#include"reaction.hpp"
#include"observable.hpp"

static double sqrt3 = 1.732050807568877;

using namespace std;

struct filenames
{
        char input_file[40];
        char temp_file[40];
        char current_temp_file[40];
        char traject_file[40];
        char bnd_tra_file[40];
        char bnd_dst_file[40];
        char average_file[40];
        char density_file[40];
        char bcd_hbk_file[40];
        char ranvar_file[40];
        char reactions_file[40];
        char observables_file[40];
};

struct parameter
{
    unsigned int what_to_do;
    unsigned int max_outputs;
    unsigned int outputs;
    unsigned int init_random_seed;
    unsigned int random_seed;
    MTRand::uint32* generator_state;
    unsigned int generator_state_len;

    double rlx_time;
    double meas_time;
    double main_time;
    unsigned long rlx_steps;
    unsigned long meas_steps;
    unsigned long main_steps;
    
    double last_time;
    double next_time;
    unsigned long ac_run_step;
    double ac_run_time;
    unsigned long ac_meas_step;
    double bgn_meas;
    double end_meas;

    unsigned int no_trajec;
    unsigned int no_inst_outputs;

    unsigned int max_rows;
    unsigned int max_cols;

    double subv_size;

    unsigned int no_nucl;
    unsigned int no_z_nucl;
    unsigned int no_r_nucl;
    unsigned int no_reac;
    unsigned int no_para;
    unsigned int no_spec;
    unsigned int no_obs;
    unsigned int no_neig;

    double* threshold;
    observable* obs;

    reaction_class* reaction;
    double* reac_parameter;
    double* diff_rate_const;
    double* configuration;

    random_var** obs_config;
    random_var** avg_obs_config;
    histo1D* obs_boundary_pos;
    histo1D* obs_boundary_no;
    
    int tree_size;
    double* tree_value;
    int* tree_index;
    double* event_time;
    
    bool output_all;
    int instant_outputs;
    
    char dump[200];

};

class run_parameter : public parameter, public filenames
{
    public:
        run_parameter(void){};
        int read_filenames(int argc, char *argv[]);
        int read_parameter(void);
        int read_temp_file(void);
        int open_output(ofstream& fout,const char * filename, bool header);
        int open_output(ofstream& fout,const char * filename, bool header, int prec);
        int open_output(ofstream& fout,const char * filename, ios_base::openmode, int prec);
        int write_header(ofstream& fout,const char * filename);
        int close_output(ofstream& fout);
};
#endif // _RUN_PARAMETER_HPP
