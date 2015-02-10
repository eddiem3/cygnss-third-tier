ORIG_INC_PATH = -I/usr/include
ORIG_LAB_PATH = -L/usr/lib

NVVM_INC_PATH = -I/usr/local/cuda/nvvm/include
NVVM_LIB_PATH = -L/usr/local/cuda/nvvm/lib64

CUDA_INC_PATH = -I/usr/local/cuda/include
CUDA_LIB_PATH = -L/usr/local/cuda/lib64




correlation:
	g++ cross_correlation.cpp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lafcuda -lnvvm -o cross_correlation.o -g $(NVVM_INC_PATH) $(NVVM_LIB_PATH) $(CUDA_INC_PATH) $(CUDA_LIB_PATH) 

retrevial: 
	pgfortran -Mcuda -o retrevial retrevial.cuf -L/usr/local/cuda/lib64 -lcufft -Mallocatable=03 -g

full:
	g++ -c cross_correlation.cpp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lafcuda -lnvvm -o cross_correlation.o -g $(NVVM_INC_PATH) $(NVVM_LIB_PATH) $(CUDA_INC_PATH) $(CUDA_LIB_PATH) 

	gfortran cross_correlation.o retrevial.f90 -o matched_filter -lafcuda


clean:
	rm -rf *.o retrevial cross_correlation.o