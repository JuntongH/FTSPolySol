NVCC ?= nvcc
TARGET ?= main
NVCCFLAGS ?= -O3 -std=c++17
CUDA_ARCH ?=
SOURCES := $(wildcard *.cu)

$(TARGET): $(SOURCES) $(wildcard *.h)
	$(NVCC) $(NVCCFLAGS) $(CUDA_ARCH) -o $@ $(SOURCES) -lcufft -lcurand

clean:
	rm -f $(TARGET) *.o
