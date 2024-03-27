//  ~~~~~~~~~~~~~~~~~~~~~~
//  analyse_trajectory.cpp
//  ~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~
//         BY TOMEK    
//         Dec 2010  
// ~~~~~~~~~~~~~~~~~~~~~~~

#include"analyse_trajectory.hpp"

#define NO_COLS 26
#define T_MIN   0.0
#define T_SHORT 6E+09       // 1E+5 for standard toggle switch, 6E+9 for rescaled
#define MAX_T 600000
#define MAX_T_CORR 5000

#define sign(x) (int)(( x > 0 ) - ( x < 0 ))


int indX(double n, double thresh_x)
{
  return (int)( n >= thresh_x );
}

int indY(double n, double thresh_y)
{
  return (int)( n <= thresh_y );
}




using namespace std;

int main(int argc, char *argv[])
{
 
  MTRand mtrand(23);
  
  // Arguments
  double t_short = atof(argv[2]);       // The maximal time of the short trajectory
  double timebinsize = atof(argv[3]);   // Binsize for the switching times histogram
  
  // The array storing the input
  double col[NO_COLS]        = {0.0};
  double last_col[NO_COLS]   = {0.0};
  
  double totalX, totalY, difXY;
  double null_time, last_time, this_time;
  double weight;
  
  bool start = true;
  
  // Histograms
  histo1D his1Ddif, his1Ddia, his1DA, his1DB, hisTimes;
  histo2D his2D;
  //double his2D2[300][300];
  double weight_sum = 0.0;
  double binsize = 5.0;
  int range1D = 500;
  int range2D = 250;
  int rangeTimes = 100;
  
  // Correlation functions
  double cX[MAX_T] = {0.0};
  double cY[MAX_T] = {0.0};  
  double dt[MAX_T] = {0.0};
  double cCorrX[MAX_T_CORR] = {0.0};
  double cCorrY[MAX_T_CORR] = {0.0};
  double cXavg, cYavg, cor_weight = 0.0, avg_weight = 0.0;
  
  double thresh_x = atof(argv[4]);      // copy number threshold to check for being in one of the exclusive switch states
  double thresh_y = -1.0*atof(argv[4]);
  
  double this_switch_time = 0.0, last_switch_time = 0.0, lst_cum = 0.0;
  double this_high_X = 0.0, this_high_Y = 0.0;
  double last_high_X = 0.0, last_high_Y = 0.0;
  
  // Output files
  ofstream his1Dfile, his2Dfile, hisTimesfile, corfile, shortfile;
  ifstream infile;
  char infilename[500], hisfilename[500], buf[500];
  
  // Some counters
  int i,j,k,t,T;
  unsigned long l = 0;
      
  // Open Infile
  strcpy(infilename, argv[1]);
  infile.open(infilename);
  if(!infile){
    
      cerr << "ERROR! Could not open infile!" << endl;
      exit(EXIT_FAILURE);
  }
  
  // Initialize histograms
  // First generate filename containing scaling factors
  // 1D histograms
  strcpy(hisfilename, infilename);
  strcat(hisfilename, ".his1D");
  his1Dfile.open(hisfilename);
  his1DA.init(0, range1D, binsize);
  his1DB.init(0, range1D, binsize);
  his1Ddif.init(0, range1D, binsize);
  his1Ddia.init(0, range1D, binsize);

  // 2D histograms
  strcpy(hisfilename, infilename);
  strcat(hisfilename, ".his2D");
  his2Dfile.open(hisfilename);
  his2D.init(0, range2D, 0, range2D, binsize);
  weight_sum = 0.0;
  
  // Flipping times histogram
  strcpy(hisfilename, infilename);
  strcat(hisfilename, ".hisT");
  hisTimesfile.open(hisfilename);
  hisTimes.init(0, rangeTimes, timebinsize);
  
  // Time correlation functions
  strcpy(hisfilename, infilename);
  strcat(hisfilename, ".cor");
  corfile.open(hisfilename);
  
  // Short trajectory file
  strcpy(hisfilename, infilename);
  strcat(hisfilename, ".short");
  shortfile.open(hisfilename);
  

  last_time = 0.0;
  this_time = 0.0;
  last_switch_time = 0.0;
  this_switch_time = 0.0;
  start = true;
  T=0;
  cout << "Reading trajectory file ..." << endl;
  while(!infile.eof())
  {
    
        // Read a new data line from the trajectory file
        for(int i=0; i<NO_COLS; i++)
        {  
          // Remember values read in the last step
          last_col[i] = col[i];
          // Read new values from infile
          infile >> col[i];

        }
        
        last_time = last_col[13];
        this_time = col[13];
        
        // Copy it to the shortened file if time<t_short
        if(this_time - null_time <= t_short){
          
            for(int i=0; i<NO_COLS; i++)    shortfile << col[i] << " ";
            shortfile << last_high_X - last_high_Y << " " << last_switch_time << " " << this_switch_time << " " << lst_cum;
            shortfile << endl;
        }
        
        // Calculate total numbers of Hb and Rep (order parameters)
        // for the last line read
        /*
        totalX = last_col[15] + 2.0*(last_col[16]+last_col[24]);
        totalY = last_col[21] + 2.0*(last_col[22]+last_col[18]);
        */
        totalX = last_col[19];
        totalY = last_col[25];
        difXY = totalX - totalY;
        // Calculate the weight of this data
        weight = this_time - last_time;
        

        if(!start){
        // First line discarded because last_cn and last_time are not defined then yet
          
            // Species number and difference histograms
            his1DA.add( (int)(totalX/binsize), weight);
            his1DB.add( (int)(totalY/binsize), weight);
            his1Ddif.add( (int)(difXY/binsize+range1D/2), weight);
            // Histogram along the diagonal
            if(totalX==totalY) his1Ddia.add( (int)(totalX/binsize), weight);
                                
            // 2D plane histogram
            his2D.add((int)(totalX/binsize), (int)(totalY/binsize), weight);
            /*
            if((int)(totalX/binsize)<range2D && (int)(totalY/binsize)<range2D)
            {
                his2D2[(int)(totalX/binsize)][(int)(totalY/binsize)] += weight;
                weight_sum += weight;
            }
            */
            
            // Save data for calculation of the time correlation function
            if(T<MAX_T){    // only if the array is not yet full
              
                assert(indX(difXY, thresh_x) + indY(difXY, thresh_y)!=2);
                cX[T] = indX(difXY, thresh_x);
                cY[T] = indY(difXY, thresh_y);
                dt[T] = weight;
                T++;
            }
            
            // Switching times histogram
            //
            // Check for switching event
            // If it just switched to the high-X state
            if( last_high_Y && indX(difXY, thresh_x) )
            {
                last_switch_time = this_switch_time;
                this_switch_time = 0.5*(this_time + last_high_Y);
                last_high_Y = 0.0;
                assert(last_high_X == 0.0);
                hisTimes.add((int)((this_switch_time - last_switch_time)/timebinsize), 1.0);
                lst_cum += (this_switch_time - last_switch_time);
            }
            // If it just switched to the high-Y state
            if( last_high_X && indY(difXY, thresh_y) )
            {
                last_switch_time = this_switch_time;
                this_switch_time = 0.5*(this_time + last_high_X);
                last_high_X = 0.0;
                assert(last_high_Y == 0.0);
                hisTimes.add((int)((this_switch_time - last_switch_time)/timebinsize), 1.0);
                lst_cum += (this_switch_time - last_switch_time);
            }
              
            // Continously track time of last high X or Y states
            if(indX(difXY, thresh_x))     last_high_X = this_time;
            if(indY(difXY, thresh_y))     last_high_Y = this_time;
            
        }
        else    null_time = this_time;
        
        // Set startflag to false after first line is read
        start = false;
        l++;
  
  } // end while(!infile.eof())
  
  cout << "Read " << l << " lines." << endl;
  
  cout << "Calculating time correlation function ..." << endl;
  if(T<MAX_T_CORR) cout << "WARNING! Not enough datapoints to calculate correlation function!" << endl;
  for(int t=0; t<=T; t++){
    
    if( t%(MAX_T/10)==0 ){
      
      cout << ".";
      cout.flush();
      
    }
      
    if(T-t > MAX_T_CORR){
          
            for(int tc=0; tc<MAX_T_CORR; tc++)
            {
                cCorrX[tc] += cX[t] * cY[t+tc];
                cCorrY[tc] += cY[t] * cX[t+tc];
            }
            
            cXavg += cX[t];
            cYavg += cY[t];
                
            cor_weight++;
        }    
  }
  cout << endl;
  
  
  cout << "Writing the output ..." << endl;
  // Output the histograms      
  // 2D histogram
  his2D.normalise();
  his2Dfile.precision(5);
  his2Dfile << "# sum of probabilities = " << his2D.sum_probab() << endl;
  for(int i=0; i<his2D.no_cols(); i++){
    
    for(int j=0; j<his2D.no_rows(); j++){
      his2Dfile << i*binsize << " "
                << j*binsize << " "
                << his2D.probab_at(i,j) << " "
                //<< his2D2[i][j]/weight_sum
                << endl;
    }
    his2Dfile << endl;
    
  }
  // 1D histograms
  his1DA.normalise();
  his1DB.normalise();
  his1Ddif.normalise();
  his1Ddia.normalise();
  his1Dfile.precision(5);
  his1Dfile << "# sum of probabilities his1DA = " << his1DA.sum_probab() << endl;
  his1Dfile << "# sum of probabilities his1DB = " << his1DB.sum_probab() << endl;
  his1Dfile << "# sum of probabilities his1Ddif = " << his1Ddif.sum_probab() << endl;
  his1Dfile << "# sum of probabilities his1Ddia = " << his1Ddia.sum_probab() << endl;
  for(int i=0; i<his1Ddif.no_bins(); i++){
      
      his1Dfile << i*binsize  << " " 
                  << his1DA.probab_at(i) << " "
                  << his1DB.probab_at(i) << " "
                  << his1Ddif.probab_at(i) << " "
                  << his1Ddia.probab_at(i) << " ";
                  
      // Size of 2D histogram might differ from 1D histogram size
      if( i<his2D.no_rows() && i<his2D.no_cols() )
                  his1Dfile << his2D.probab_at(i,i) << endl;
      else
                  his1Dfile << endl;
  }
  // Switching time histogram
  hisTimes.normalise();
  hisTimesfile.precision(5);
  hisTimesfile << "# sum of probabilities hisTimes  = " << hisTimes.sum_probab() << endl;
  for(int i=0; i<hisTimes.no_bins(); i++){
      
      hisTimesfile << i*timebinsize  << " " 
                  << hisTimes.probab_at(i)  << " "
                  << (hisTimes.probab_at(i)>0.0 ? log(hisTimes.probab_at(i)) : 0.0) << " "
                  << endl;
  }
  
  // Correlation times file
  for(int tc=0; tc<MAX_T_CORR; tc++)
       corfile << tc << " "
               << cCorrX[tc]/cor_weight << " "
               << cXavg/cor_weight << " "
               << cCorrX[tc]/cXavg << " "
               << cCorrX[tc] << " "
               << cCorrY[tc] << " "
               << cXavg << " "
               << cYavg << " "
               << cor_weight
               << endl;
  
  // Clean up
  his1Dfile.flush();
  his1Dfile.close();
  his2Dfile.flush();
  his2Dfile.close();
  hisTimesfile.flush();
  hisTimesfile.close();
  corfile.flush();
  corfile.close();
  
  cout << "Done." << endl << endl;

  
} // end main
