#include <iostream>
#include <chrono>
#include <cmath>
#include <omp.h>

using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
int Nth; 
unsigned int count = 0;

int minthreads, maxthreads;

void myprimes(const unsigned int tot) {
   double start_time, run_time;
   milliseconds tot_ms;
   if(tot > 1) count++;
   for (int nthread = minthreads ; nthread <= maxthreads ; nthread++){
      omp_set_num_threads(nthread);
      auto t0 = high_resolution_clock::now();
      start_time = omp_get_wtime();
      count = 0;
      #pragma omp parallel for reduction(+:count)
      for(unsigned int i=3; i<=tot; i++){
         Nth = omp_get_num_threads();
         unsigned int s=sqrt(i);
         bool isprime=true;
         if(i%2 == 0) isprime=false;
         for(unsigned int j=3; (j<=s) && (isprime); j+=2) if (i%j == 0) isprime=false;
         if(isprime) {
            count++;
         }
      }
      auto t1 = high_resolution_clock::now();
      milliseconds tot_ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
      run_time = omp_get_wtime() - start_time;
      tot_ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
      std::cout << "##########################################################\n";
      std::cout << "time OMP in seconds is " << run_time << std::endl;
      //std::cout << "high_resolution_clock milliseconds = " << tot_ms.count() << std::endl;
      std::cout << "N Threads = " << Nth << std::endl;
      std::cout << "N Prime  = " << count << std::endl;

   }
}

int main(int argc, char* argv[])
{
 if (argc < 4) {
  std::cerr << "Usage: " << argv[0] << " <MIN THREADS> <MAX THREADS> <NUMBER>" << std::endl;
  return 1;
 }
 minthreads = atoi(argv[1]);
 maxthreads = atoi(argv[2]);
 unsigned int tot = atoi(argv[3]);
 myprimes(tot);
 //milliseconds tot_ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
 //std::cout << "milliseconds = " << tot_ms.count() << std::endl;
 //std::cout << "N Threads = " << Nth << std::endl;
 //std::cout << "N Prime  = " << count << std::endl;
return 0;
}
