ORIG_INC_PATH = -I/usr/include
ORIG_LAB_PATH = -L/usr/lib

CUFFTW_INC_PATH = -I/usr/local/cuda-6.0/targets/x86_64-linux/include
CUFFTW_LIB_PATH = -L/usr/local/cuda-6.0/targets/x86_64-linux/lib

all: fftw_orig fftw_cuda

fftw_orig: fftw_example.c
	gcc -g -o $@ fftw_example.c -lfftw3 -lfftw3f -lm $(ORIG_INC_PATH) $(ORIG_LIB_PATH)

fftw_cuda: fftw_example.c
	gcc -g -o $@ fftw_example.c -DCUFFTW -lcufftw -lcufft -lm $(CUFFT_INC_PATH) $(CUFFT_LIB_PATH)

retrevial: retrevialfft.cu
	nvcc -o $@ retrevialfft.cu -lcuda -lcudart -lcufft -lm $(CUFFT_INC_PATH) $(CUFFT_LIB_PATH)

retrevial_fortran: 
	pgfortran -Mcuda -o retrevial retrevial.cuf -L/usr/local/cuda/lib64 -lcufft -Mallocatable=03 -g

clean:
	rm -rf *.o fftw_orig fftw_cuda retrevial retrevial_fortran