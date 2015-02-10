#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <ctime>

//Array Fire Includes
#include <arrayfire.h>
#include <af/util.h>

using namespace af;

/* A function to compute the cross correlation of two signals.
   Useful for finding the correlation of a waveform and a storm signal
   @param s1 - The first signal
   @param storm - The second signal
*/
extern "C" void correlateSignals(double s1[], double s2[])
{

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

  array result = numerator/denominator;

  print("final", result);
  
  //print("signal", signal);
    
    
  // //return result;

 
}

//int main(int argc, char ** argv)
//{
// //Read in Nature Run Data
// std::cout << "Reading interpolated nature run data" << "\n";
// std::ifstream nature("interpolatedNatureRun.txt");

// std::vector<array> nature_vectors;

// long double n;
// while(!nature.eof())
//   {
//     array nature_signal(600,1);

//     for(int i = 0; i < nature_signal.dims(0); i++)
// {
//   nature >> n;
//   nature_signal(i) = n;  
//   //std::cout << n << "\n";
// }
//     nature_vectors.push_back(nature_signal);
//   }

// std::cout << nature_vectors.size() << "\n";

// print("a",nature_vectors.at(1));


//Read in Model Run Data
// std::cout << "Reading interpolated model run data" << "\n";
// std::ifstream model("Willoughby56_530-w.bin2.bin", std::ios::binary | std::ios::in);

// std::vector<array> model_vectors;

//  array model_signal(600,1);
//  long double o;
//  double duration;

//  std::clock_t start;

// // start = std::clock();
// while(!model.eof())
//   {

//     for(int i = 0; i < model_signal.dims(0); i++)
// {
//   model >> o;
//   model_signal(i) = o;  
//   //std::cout << n << "\n";
// }
//     model_vectors.push_back(model_signal);
//   }
// duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

// //std::cout << model_vectors.size() << "\n";

// // print("a",model_vectors.at(1));

// std::cout<<"Subarrays:  "<< duration <<'\n';

// std::cout << model_vectors.size() << "\n";


//  *getWilloughby_();



// try {
    
//   //Select a device and display info
//   int device = argc > 1 ? atoi(argv[1]) : 0;
//   deviceset(device);
//   info();
    
//   timer::start();
//   for(int i=0; i < nature_vectors.size(); i++)
//     {
// for(int j=0; j < 55661; j++)
//   {
//     correlateSignals(nature_vectors[i], model_vectors[j]);
//   }
//     }
//   timer::stop();

    
// } catch (af::exception& e) {
//   fprintf(stderr, "%s\n", e.what());
// }  
//  return 0;
//}
