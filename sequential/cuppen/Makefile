VERSION = 4
EXECNAME = cuppen_v$(VERSION)

SRCS := cuppen_v4.c daq_v3.c merge_v2.c deflate_v4.c rank_one_v2.c zerodandc_n_v2.c calculate_vi_v2.c calculate_normq_v2.c calculate_q_v2.c mult_v1.c mag_v4.c gen_v2.c err_v3.c hib_time_v2.c io_v2.c
OBJS := $(SRCS:.c=.o)

CU_SRCS :=
CU_OBJS := $(CU_SRCS:.cu=.o)

CUDADIR ?= /opt/cuda
MAGMADIR = /home/equipo2/maestria/magma-1.5.0-beta2

MAGMA_CFLAGS := -DADD_ -DHAVE_CUBLAS $(shell pkg-config --cflags magma)
MAGMA_LIBS :=  $(shell pkg-config --libs magma)

CC = /usr/lib64/openmpi/1.4-gcc/bin/mpicc

DEFS =

INC_PATHS =
LIB_PATHS =

CFLAGS = -O3 -m64 -fpic -fopenmp -Wall
CPPFLAGS = $(DEFS) $(INC_PATHS)
LDFLAGS = $(LIB_PATHS)
LDLIBS = -lgomp -lrt -lc

NVCC = nvcc

NVDEFS =

NVINC_PATHS = -I$(CUDADIR)/include
NVLIB_PATHS = -L$(CUDADIR)/lib64

NVCFLAGS = -O3 -m64 -prec-sqrt=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20 -Xcompiler "-Wall -O3 -m64 -fopenmp -fPIC"
NVCPPFLAGS = $(NVDEFS) $(NVINC_PATHS)
NVLDFLAGS = $(NVLIB_PATHS)
NVLDLIBS = -lcudart -lcublas -lm -lc

.PHONY: all clean cleanobjs

.DEFAULT_GOAL := all

all : $(EXECNAME)
#	$(RM) $(OBJS) $(CU_OBJS)

$(EXECNAME) : $(OBJS) $(CU_OBJS); @echo 
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS) $(MAGMA_LIBS) -o $(EXECNAME)
	@echo 

$(OBJS) : $(SRCS)

$(CU_OBJS) : $(CU_SRCS)

%.o: %.c; @echo 
	$(CC) $(CPPFLAGS) $< $(CFLAGS) $(MAGMA_CFLAGS) -c -o $@

%.o: %.cu; @echo 
	$(NVCC) $(NVCPPFLAGS) $< $(NVCFLAGS) $(NVLDFLAGS) $(NVLDLIBS) -c -o $@

cleanobjs:
	$(RM) $(OBJS) $(CU_OBJS)

clean :
	$(RM) $(OBJS) $(CU_OBJS) $(EXECNAME)

