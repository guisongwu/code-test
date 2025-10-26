# compiler settings 
PETSC_DIR=/opt/petsc-3.23.2
PETSC_ARCH=
MPI_DIR=/usr/lib/x86_64-linux-gnu/openmpi
EIGEN_DIR=/usr/include/eigen3

# include ${PETSC_DIR}/lib/petsc/conf/petscvariables
# include ${PETSC_DIR}/lib/petsc/conf/variables
# include ${PETSC_DIR}/lib/petsc/conf/rules

# INCLUDE_DIR_ARCH=/usr/include/x86_64-linux-gnu 
# # gmp.h, cblas.h
# LIB_DIR_ARCH=/usr/lib/x86_64-linux-gnu
# # libgmp.so, libmpfr.so, libblas.so, liblapack.so
# INCLUDE_DIR=/usr/include
# # mpfr.h, mpreal.h, lapack.h, lapacke.h
# LIB_DIR=/usr/lib

CPROGRAM_DIR=/home/wugs/cprogram

# CC = mpicc 
# CXX = mpicxx 
CC = gcc 
CXX = g++ 

C_SRCS = $(wildcard *.c)
CPP_SRCS = $(wildcard *.cpp)
OBJS = $(C_SRCS:.c=.o) $(CPP_SRCS:.cpp=.o)

CFLAGS   = -g -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -std=gnu17
CXXFLAGS = -g -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -Wno-unused-but-set-variable -Wno-comment -Wno-sign-compare -std=c++20

# C/C++ PreProcessor options
CPPFLAGS = -I/usr/include/eigen3 -I$(CPROGRAM_DIR) -I$(PETSC_DIR)/include -I$(MPI_DIR)/include
# Library searching dir options
LDFLAGS = -L$(MPI_DIR)/lib -Wl,-rpath,$(MPI_DIR)/lib -L$(PETSC_DIR)/lib -Wl,-rpath,$(PETSC_DIR)/lib
# Library link options
LDLIBS = -lm -lpetsc -lmpi -llapacke -llapack -lblas


# target files 
TARGETS = testcpp poisson-rt0 helmholtz-rt0 helmholtz-rt0-sphere helmholtz-rt1 helmholtz-rt1-sphere mesh poisson-rt1 poisson-rt0-mix poisson-rt0-sphere poisson-rt1-sphere
CTARGETS = poisson testc testpetsc

# default target(when "make" is called, all the target files will be built.)
all: $(TARGETS)

########################### Link ################################
# Generic link rule for c program single .o -> executable
$(CTARGETS): %: %.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)
# Generic link rule for cpp program single .o -> executable
$(TARGETS): %: %.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)


########################### Compile #############################
# C program compile rules
$(filter %.o,$(C_SRCS:.c=.o)): %.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

# C++ program compile rules
$(filter %.o,$(CPP_SRCS:.cpp=.o)): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@



# general build rule for C program to object file
# %.o: %.c
# 	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@
# general build rule for CPP program to object file
# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@



# clean object and target files
clean:
	rm -f *.o $(TARGETS) $(CTARGETS)

# 伪目标声明
.PHONY: all clean
