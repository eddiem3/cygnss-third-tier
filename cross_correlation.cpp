#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

//Array Fire Includes
#include <arrayfire.h>
#include <af/util.h>

using namespace af;


/* A function to compute the cross correlation of two signals.
   Useful for finding the correlation of a waveform and a storm signal
   @param signal - The waveform database array
   @param storm - The real storm array 
*/
static array correlateSignals(array signal, array storm)
{

  array fourierSignal;
  array fourierStorm;

  //compute the fft of the signal in the waveform database
  fourierSignal = fft(signal);  

  //print("fdmims", fourierSignal.dims());

  //compute the complex conjugate of the fft of the storm 
  fourierStorm = conjg(fft(storm));
  

  //make sure that the variables associated with each operation are evaulated
  fourierSignal.eval();
  fourierStorm.eval();

  
  //magnitude 
  //sum and multiply signals
  array numerator = conjg(sum(fourierSignal * fourierStorm)) * sum(fourierSignal * fourierStorm);

  array denominator = sum(conjg(fourierSignal) * fourierSignal) * sum(conjg(fourierStorm) * fourierStorm);

  array result = numerator/denominator;

  print("final", result);

    
  return result;


}

int main(int argc, char ** argv)
{
  try {

    //Select a device and display info
    int device = argc > 1 ? atoi(argv[1]) : 0;
    deviceset(device);
    info();


    //    timer start = timer::start();

    //for(int i=1; i < 600; i++)
    // {
	array A = randu(600,1, f32);
	array B = randu(600,1, f32);
	correlateSignals(A,A);
	//}
	//std::cout << timer::stop(start) << std::endl;
            
  } catch (af::exception& e) {
    fprintf(stderr, "%s\n", e.what());
  }
  return 0;
}
