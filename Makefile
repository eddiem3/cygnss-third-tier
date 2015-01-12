ORIG_INC_PATH = -I/usr/include
ORIG_LAB_PATH = -L/usr/lib

#CUFFTW_INC_PATH = -I/usr/local/cuda-6.0/targets/x86_64-linux/include
#CUFFTW_LIB_PATH = -L/usr/local/cuda-6.0/targets/x86_64-linux/lib


correlation:
	g++ cross_correlation.cpp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lafcuda -o cross_correlation.o -g

retrevial: 
	pgfortran -Mcuda -o retrevial retrevial.cuf -L/usr/local/cuda/lib64 -lcufft -Mallocatable=03 -g

clean:
	rm -rf *.o retrevial 