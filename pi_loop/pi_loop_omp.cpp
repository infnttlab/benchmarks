/*

This program will numerically compute the integral of

                  4/(1+x*x) 
				  
from 0 to 1.  The value of this integral is pi -- which 
is great since it gives us an easy way to check the answer.

History: Original C version written by Tim Mattson, 11/99.
         C++ version by A.Ferraro, D.Cesini 5/2014
*/
#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <chrono>

using namespace std;

//static long num_steps = 2000000000;
double step;

int main (int argc, char* argv[])
{
 if (argc < 4) {
      cerr << "Usage: " << argv[0] << " <THREADS> <STEPS> <VERBOSE>" << std::endl;
      return 1;
  }
 int max_threads = atoi(argv[1]);
 long num_steps = atol(argv[2]);
 int verbose = atoi(argv[3]);
 int n_threads;
 double x, pi, sum = 0.0;
 double start_time, run_time;
 chrono::high_resolution_clock::time_point t01, t02;
 chrono::duration<double> dt1;
 step = 1.0/(double) num_steps;
 cout << "PI calculator (" << argv[0] << ")" << endl;
 for (n_threads = 1; n_threads<=max_threads; n_threads++){
  sum = 0.0;
  omp_set_num_threads(n_threads);
  start_time = omp_get_wtime();
  t01 = chrono::high_resolution_clock::now();
  #pragma omp parallel
  {
   double x;
   #pragma omp single
   {if(verbose)
     cout << "num_threads: " << omp_get_num_threads() << "  ";
   }
   #pragma omp for reduction(+:sum)
   for (int i=1;i<= num_steps; i++){
     x = (i-0.5)*step;
     sum = sum + 4.0/(1.0+x*x);
    }
   }
   pi = step * sum;
   run_time = omp_get_wtime() - start_time;
   t02 = chrono::high_resolution_clock::now();
   dt1 = chrono::duration_cast<chrono::duration<double>>(t02 - t01);
   if(verbose)
    cout << "PI: " << pi << "  Time:  " << dt1.count() <<  " secs" << endl; 
   else
    cout << dt1.count() << endl;
 }
}	  

