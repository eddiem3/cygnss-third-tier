#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

//Array Fire Includes
#include <arrayfire.h>
#include <af/util.h>

using namespace af;

//extern {
// correlateSignals_()
//}


/* A function to compute the cross correlation of two signals.
   Useful for finding the correlation of a waveform and a storm signal
   @param signal - The waveform database array
   @param storm - The real storm array 
*/
static array correlateSignals(float signalArray[], float stormArray[])
{


  //Convert the storm signals into Array Fire array type
  array fourierSignal(600,1,signalArray);
  array fourierStorm(600,1,stormArray);
 
  //compute the fft of the signal in the waveform database
  fourierSignal = fft(fourierSignal);  

  //print("fdmims", fourierSignal.dims());

  //compute the complex conjugate of the fft of the storm 
  fourierStorm = conjg(fft(fourierStorm));
  

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
    
    float signalA[600];
    float signalB[600];

    for(int i = 0; i < 600; i++)
      {
	signalA[i] = i;
	signalB[i] = i;
      }

    // A = randu(600,1, f32);
    //array B = randu(600,1, f32);
    correlateSignals(signalA,signalB);
    
  }             catch (af::exception& e) {
    fprintf(stderr, "%s\n", e.what());
  }
  return 0;
}
