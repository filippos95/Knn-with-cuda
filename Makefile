CMAIN = ex3

parallel: nvcc -arch=sm_70 -o cuda-knn-parallel main-cuda.cu -run



clean:
	rm -f cuda-knn-parallel cuda-knn-serial
