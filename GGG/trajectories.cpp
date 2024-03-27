// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// t r a j e c t o r i e s . c p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "trajectories.hpp"

using namespace std;

int trajectory(run_parameter input)
{
	// random number generator
	MTRand::uint32 bigSeed[4] = { 0x123, 0x234, 0x345, 0x456 };
	MTRand mtrand( bigSeed , 4 );

	// open output file
	ofstream trac; input.open_output(trac,input.traject_file,0);

	// create system
	cortex fly; fly.init_cortex(mtrand,input);

	// relaxation run
	unsigned long ac_rlx_step = 0;
	double ac_rlx_time = 0.0;
	do {
		if(ac_rlx_step >= input.rlx_steps) break;
		ac_rlx_time = fly.nuclear_queue.get_value(0);
		fly.gillespie_step(mtrand);
		ac_rlx_step++;
	} while(ac_rlx_time < input.rlx_time);

	// production run
	unsigned long ac_max_step = 0;
	double last_time = ac_rlx_time;
	double next_time = fly.nuclear_queue.get_value(0);
	if(input.meas_time == 0.0 || input.meas_steps == 0) // choose : no averaging
	{
		do {
			for(unsigned int r_nucl = 0; r_nucl < input.no_r_nucl; r_nucl++)
			{
				double r_pos = double(r_nucl+0.5)/double(input.no_r_nucl);
				double r_dst = double(r_nucl+0.5)*double(input.subv_size);
				for(unsigned int z_nucl = 0; z_nucl < input.no_z_nucl; z_nucl++)
				{
					double z_pos = double(z_nucl+0.5)/double(input.no_z_nucl);
					double z_dst = double(z_nucl+0.5)*double(input.subv_size);
					trac << input.subv_size << " ";
					trac << input.no_r_nucl << " " << r_nucl << " " << r_pos << " " << r_dst << " ";
					trac << input.no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " ";
					trac << ac_rlx_step << " " << ac_max_step << " " << ac_rlx_time << " ";
					trac << 0.5*(next_time+last_time) << " " << last_time << " ";
					trac << double(fly.get_config(r_nucl,z_nucl,0)) << " ";
					trac << double(fly.get_config(r_nucl,z_nucl,1)) << " ";
					trac << double(fly.get_config(r_nucl,z_nucl,7)) << " ";
					trac << endl;
				}
			}
			fly.gillespie_step(mtrand);

			last_time = next_time;
			next_time = fly.nuclear_queue.get_value(0);

			ac_max_step++;
			if (ac_max_step >= input.main_steps) break;

		} while(next_time <= input.main_time+ac_rlx_time);
	}
	else // choose : average
	{
		do {
			// temporary storage
			double * tmp_bcd_config = new double[input.no_r_nucl*input.no_z_nucl];
			double * tmp_hbk_config = new double[input.no_r_nucl*input.no_z_nucl];
			double * tmp_obf_config = new double[input.no_r_nucl*input.no_z_nucl];
			for(unsigned int r_nucl = 0; r_nucl < input.no_r_nucl; r_nucl++)
			{
				for(unsigned int z_nucl = 0; z_nucl < input.no_z_nucl; z_nucl++)
				{
					tmp_bcd_config[z_nucl+r_nucl*input.no_z_nucl] = 0.0;
					tmp_hbk_config[z_nucl+r_nucl*input.no_z_nucl] = 0.0;
					tmp_obf_config[z_nucl+r_nucl*input.no_z_nucl] = 0.0;
				}
			}
			// begin averaging
			unsigned long ac_avg_step = 0;
			double bgn_avg = last_time;
			double end_avg = last_time + input.meas_time;
			do {
				// add configuration to temporary storage
				for(unsigned int r_nucl = 0; r_nucl < input.no_r_nucl; r_nucl++)
				{
					for(unsigned int z_nucl = 0; z_nucl < input.no_z_nucl; z_nucl++)
					{
						tmp_bcd_config[z_nucl+r_nucl*input.no_z_nucl] += double(fly.get_config(r_nucl,z_nucl,0))*(next_time-last_time);
						tmp_hbk_config[z_nucl+r_nucl*input.no_z_nucl] += double(fly.get_config(r_nucl,z_nucl,1))*(next_time-last_time);
						tmp_obf_config[z_nucl+r_nucl*input.no_z_nucl] += double(fly.get_config(r_nucl,z_nucl,7))*(next_time-last_time);
					}
				}
				fly.gillespie_step(mtrand);

				last_time = next_time;
				next_time = fly.nuclear_queue.get_value(0);

				ac_avg_step++;
				if (ac_avg_step >= input.meas_steps) break;

			} while (last_time < end_avg); // end averaging

			ac_max_step += ac_avg_step;

			// print averaged configuration
			for(unsigned int r_nucl = 0; r_nucl < input.no_r_nucl; r_nucl++)
			{
				double r_pos = double(r_nucl+0.5)/double(input.no_r_nucl);
				double r_dst = double(r_nucl+0.5)*double(input.subv_size);
				for(unsigned int z_nucl = 0; z_nucl < input.no_z_nucl; z_nucl++)
				{
					double z_pos = double(z_nucl+0.5)/double(input.no_z_nucl);
					double z_dst = double(z_nucl+0.5)*double(input.subv_size);
					trac << input.subv_size << " ";
					trac << input.no_r_nucl << " " << r_nucl << " " << r_pos << " " << r_dst << " ";
					trac << input.no_z_nucl << " " << z_nucl << " " << z_pos << " " << z_dst << " ";
					trac << ac_rlx_step << " " << ac_max_step << " " << ac_rlx_time << " ";
					trac << 0.5*(next_time+last_time) << " " << last_time << " ";
					trac << tmp_bcd_config[z_nucl+r_nucl*input.no_z_nucl]/double(last_time-bgn_avg) << " ";
					trac << tmp_hbk_config[z_nucl+r_nucl*input.no_z_nucl]/double(last_time-bgn_avg) << " ";
					trac << tmp_obf_config[z_nucl+r_nucl*input.no_z_nucl]/double(last_time-bgn_avg) << " ";
					trac << endl;
				}
			}

			delete [] tmp_bcd_config;
			delete [] tmp_hbk_config;
			delete [] tmp_obf_config;

			if (ac_max_step >= input.main_steps) break;

		} while(next_time <= input.main_time+ac_rlx_time);
	} // end if : averaging
	input.close_output(trac);

	return(EXIT_SUCCESS);
}

