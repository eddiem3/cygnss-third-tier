#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <ctime>
#include <typeinfo>

//Array Fire Includes
#include <arrayfire.h>
#include <af/util.h>

using namespace af;

/* A function to compute the cross correlation of two signals.
   Useful for finding the correlation of a waveform and a storm signal
   @param s1 - The first signal
   @param storm - The second signal
*/
extern "C" double correlateSignals(double s1[], double s2[])
{


  try {
    
    //Select a device and display info
    //int device = argc > 1 ? atoi(argv[1]) : 0;
    //deviceset(0);


    array signal(600,1,s1);
    array storm(600,1, s2);
  
    //compute the fft of the signal in the waveform database
    array fourierSignal = fft(signal);  

    //compute the complex conjugate of the fft of the storm 
    array fourierStorm = conjg(fft(storm));
  
    //make sure that the variables associated with each operation are evaulated
    fourierSignal.eval();
    fourierStorm.eval();

    //magnitude 
    //sum and multiply signals
    array numerator = conjg(sum(fourierSignal * fourierStorm)) * sum(fourierSignal * fourierStorm);

    array denominator = sum(conjg(fourierSignal) * fourierSignal) * sum(conjg(fourierStorm) * fourierStorm);

    //double num = numerator(0,0).scalar<double>();

    //print("num" ,real(numerator(0,0)));

    array result = numerator/denominator;
    
    double correlation = real((numerator/denominator)(0,0)).scalar<double>();




  } catch (af::exception& e) {
    fprintf(stderr, "%s\n", e.what());
  }  

  return correlation;
}
 
