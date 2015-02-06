ORIG_INC_PATH = -I/usr/include
ORIG_LAB_PATH = -L/usr/lib

CUDA_INC_PATH = -I/usr/local/cuda-6.0/targets/x86_64-linux/include
CUDA_LIB_PATH = -L/usr/local/cuda-6.0/targets/x86_64-linux/lib


correlation:
	g++ cross_correlation.cpp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lafcuda -o cross_correlation.o -g $(CUDA_INC_PATH) $(CUDA_LIB_PATH)

retrevial: 
	pgfortran -Mcuda -o retrevial retrevial.cuf -L/usr/local/cuda/lib64 -lcufft -Mallocatable=03 -g

full:
	g++ -c cross_correlation.cpp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lafcuda -o cross_correlation.o -g $(CUDA_INC_PATH) $(CUDA_LIB_PATH)

	gfortran -c load_files.f90

	g++ -o matched_filter load_files.o cross_correlation.o -L/usr/local/cuda/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -lgfortran -lcudart -lcublas -lcufft -lafcuda $(CUDA_INC_PATH) $(CUDA_LIB_PATH)


clean:
	rm -rf *.o retrevial cross_correlation.o