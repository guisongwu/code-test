# compiler settings 
CC = gcc
CXX = g++
CFLAGS = -g -Wall -Wno-unused-variable -Wno-unused-parameter -Wextra -std=gnu17
CXXFLAGS = -g -Wall -Wno-unused-variable -Wno-unused-parameter -Wextra -std=c++20
LDFLAGS = -L/usr/local/lib -Wl,-rpath=/usr/local/lib
LDLIBS = -lm

# target files 
TARGETS = main poisson learn testc testcpp

# default target(when "make" is called, all the target files will be built.)
all: $(TARGETS)

########################### Link ################################

# build rule of "testc"
testc: testc.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# build rule of "testcpp"
testcpp: testcpp.o
	$(CXX) $(CXXFLAGS) $^ -o $@

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
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean object and target files
clean:
	rm -f *.o $(TARGETS)

# 伪目标声明
.PHONY: all clean
