//  ~~~~~~~~~~~~~~~~~~~~~~~~~
//  g i l l e s p i e . c p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~

#include<cfloat>
#include<cstdlib>
#include<string.h>
#include <math.h>
#include<algorithm>
#include<iostream>
#include<fstream>

#define MAX_T 600000
#define MAX_T_CORR 10000

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2) 
    {
        cerr << "You have to pass an argument to function main!" << endl;
        exit(EXIT_FAILURE);
    }
    
    ifstream infile;
    ofstream outfile;
    
    double cA[MAX_T] = {0.0};
    double cB[MAX_T] = {0.0};  
    double t[MAX_T] = {0.0};
    double cCorrA[MAX_T_CORR] = {0.0};
    double cCorrB[MAX_T_CORR] = {0.0};
    double cAavg  = 0.0;
    double cBavg  = 0.0;
    double weight = 0.0;
    
    unsigned long long tc=0, tmax=0, T, max_t_corr = MAX_T_CORR;
    int i,j,k,l;
    
    double buf[50];
    char charbuf[110];
    
    
    if(!infile)
    {
        cerr << "Specified infile not found!" << endl;
        exit(EXIT_FAILURE);
    }
    infile.open(argv[1]);
    infile.precision(8);
    cout.precision(8);
    
    
    // Read and discard first 10 lines
    for(int i=0; i<10; i++){
      
        infile.getline(charbuf,100);
    }
    // Now process the rest of the file
    tc = 0;
    while(!infile.eof() && tc<MAX_T)
    {
      
      for(int i=0; i<20; i++)   infile >> buf[i];

      t[tc] =  buf[0];
      cA[tc] = buf[18];
      cB[tc] = buf[19];
      
      tc++;
      
    }
    tmax = tc;
    
    // Calculate the correlation function and the time average of cA
    for(T=0; T<tmax; T++){
      
        if(tmax-T > max_t_corr){
          
            for(tc=0; tc<max_t_corr; tc++)
            {
                cCorrA[tc] += cA[T] * cB[T+tc];
                cCorrB[tc] += cB[T] * cA[T+tc];
            }
            
            cAavg += cA[T];
            cBavg += cB[T];
            weight++;
        }
    }
    
    // Output it
    for(tc=0; tc<max_t_corr; tc++)
        cout << tc << " " << cCorrA[tc]/weight << " " << cAavg/weight << " " << cCorrA[tc]/cAavg << " "
             << cCorrA[tc] << " " << cCorrB[tc] << " " << cAavg << " " << cBavg << " " << weight << endl;
    
    return(EXIT_SUCCESS);
}
