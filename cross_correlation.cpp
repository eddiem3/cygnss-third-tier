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

  //compute the complex conjugate of the fft of the storm 
  fourierStorm = conjg(fft(storm));
  

  //make sure that the variables associated with each operation are evaulated
  fourierSignal.eval();
  fourierStorm.eval();

  //multiply the signals
  array crossResult = ifft(fourierSignal * fourierStorm); 

  crossResult.eval();
  
  print ("crossResult", crossResult);
  return crossResult;
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
	array A = randu(1,600, f32);
	array B = randu(1,600, f32);
	correlateSignals(A,B);
	//}
	//std::cout << timer::stop(start) << std::endl;
    



        
  } catch (af::exception& e) {
    fprintf(stderr, "%s\n", e.what());
  }
  return 0;
}
