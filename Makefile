# compiler settings 
PETSC_DIR=/opt/petsc-3.23.2
PETSC_ARCH=

EIGEN_DIR=/usr/include/eigen3
# include ${PETSC_DIR}/lib/petsc/conf/petscvariables
# include ${PETSC_DIR}/lib/petsc/conf/variables
# include ${PETSC_DIR}/lib/petsc/conf/rules

LAPACK_DIR=/opt/OpenBLAS-0.3.29
MPREAL_DIR=
MPFR_DIR=/usr/lib/x86_64-linux-gnu
GMP_DIR=/usr/lib/x86_64-linux-gnu


CC = mpicc 
CXX = mpicxx 
CFLAGS = -g -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -Wextra -std=gnu17
CXXFLAGS = -g -Wno-unused-variable -Wno-unused-parameter -Wextra -std=c++20
CPPFLAGS = -I/usr/include/eigen3 -I$(PETSC_DIR)/include -I$(MPI_DIR)/include
LDFLAGS = -L/usr/local/lib -Wl,-rpath=/usr/local/lib -L$(MPI_DIR)/lib -Wl,-rpath=$(MPI_DIR)/lib -L$(PETSC_DIR)/lib -Wl,-rpath=$(PETSC_DIR)/lib
LDLIBS = -lm -lmpi -lpetsc


# target files 
TARGETS = main poisson learn testc testcpp poisson-rt

# default target(when "make" is called, all the target files will be built.)
all: $(TARGETS)

########################### Link ################################

# build rule of "testc"
testc: testc.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# build rule of "testcpp"
testcpp: testcpp.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# build rule of "poisson-mixed-quad"
poisson-mixed-quad: poisson-mixed-quad.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@
# poisson-mixed-quad: poisson-mixed-quad.o
#     $(CLINKER) -o $@ $^ $(PETSC_LIB)

# build rule of "main" 
main: main.o utils.o my-math.o
	$(CC) $(CFLAGS) $^ -o $@

# build rule of "poisson"
poisson: poisson.o
	$(CC) $(CFLAGS) $^ -o $@

# build rule of "learn" (cpp program)
learn: learn.o
	$(CXX) $(CXXFLAGS) $^ -o $@


######################### Compile #############################

# general build rule for C program to object file
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# general build rule for CPP program to object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# clean object and target files
clean:
	rm -f *.o $(TARGETS)

# 伪目标声明
.PHONY: all clean
