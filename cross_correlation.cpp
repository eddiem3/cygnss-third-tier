#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdio.h>

//Array Fire Includes
#include <arrayfire.h>
#include <af/util.h>

using namespace af;


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
  //Read in Nature Run Data
  std::cout << "Reading interpolated nature run data" << "\n";
  std::ifstream nature("interpolatedNatureRun.txt");

  std::vector<array> nature_vectors;

  long double n;
  while(!nature.eof())
    {
      array nature_signal(600,1);
	
      for(int i = 0; i < nature_signal.dims(0); i++)
	{
	  nature >> n;
	  nature_signal(i) = n;	  
	  //std::cout << n << "\n";
	}
      nature_vectors.push_back(nature_signal);
    }

  std::cout << nature_vectors.size() << "\n";

  print("a",nature_vectors.at(1));


  //  for(int i=0; i < nature_vectors.size(); i++){
  //print("a",nature_vectors[i]);
  //  }



  // for(std::vector<int>::const_iterator it = nature_vectors.begin(); 
  //     it != nature_vectors.end(); ++it)
  //   {
  //     print("a",*it);
  //   }



  

  //Convert std::vectors into arrays
  //Note: The willoughby data is VERY LARGE
  //Allocate on the heap instead of the stack or you
  //will get a segmentation fault
  //std::cout << "Converting will to an array" << "\n";
  //float *willoughby = (float *)malloc(33996600 * sizeof(float));
  //std::copy(will.begin(), will.end(), willoughby);
  
  //std::cout << "Converting nature to an array" << "\n";
  //float natureRun[12444];
  //std::copy(bob.begin(), bob.end(),natureRun);

  
  
  


  
 




  // std::cout << "Splitting big willougby array into subarrays" << "\n";
  // for(int i =0; i < 33996600; i++)
  //   {
  //     float data[600];

  //     for(int j = 0; j < 600; j++)
  // 	{
  // 	  data[j] = willoughby[i + 600 * j];
  // 	}	   
  //     array data_array(600,1,data);
  //     willougbyModels.push_back(data_array);      
  //   }

  // std::cout << "splitting up nature array" << "\n";
  // for(int i=0; i < 12444; i++)
  //   {
  //     float data[600];

  //     for(int j=0; j < 600; j++)
  // 	{
  // 	  data[j] = natureRun[i + 600 * j];
  // 	}

  //     array data_array(600,1,data);
  //     natureRunData.push_back(data_array);

  //   }
 

  
  // for (int i=0; i<natureRunData.size(); i++) 
  //   {
  //     print("a",natureRunData[i]);
  //   } 

	  
  //Loop through each array in the storm vector 
  //Multiply it with each of the willougby
  //  for(auto &i : natureRunData)
  //{
  //  i * i;
  //}
  






  //DON'T DELETE
  //print to std out
  // std::cout << "numbers read in:\n";
  // std::copy(bob.begin(), bob.end(), 
  //           std::ostream_iterator<float>(std::cout, " "));
  // std::cout << std::endl;


  //



  // try {
    
  //   //Select a device and display info
  //   int device = argc > 1 ? atoi(argv[1]) : 0;
  //   deviceset(device);
  //   info();
    
  //   float signalA[600];
  //   float signalB[600];

  //   for(int i = 0; i < 600; i++)
  //     {
  // 	signalA[i] = i;
  // 	signalB[i] = i;
  //     }

  //   // A = randu(600,1, f32);
  //   //array B = randu(600,1, f32);
  //   correlateSignals(signalA,signalB);
    
  // }             catch (af::exception& e) {
  //   fprintf(stderr, "%s\n", e.what());
  // }

  //  free(willoughby);
  return 0;
}
