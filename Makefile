# 编译器设置
CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -std=c11
CXXFLAGS = -Wall -Wextra -std=c++11

# 目标文件
TARGETS = main poisson-fem learn

# 默认目标（执行 make 时默认构建所有目标）
all: $(TARGETS)

# main 的构建规则
main: main.o utils.o math.o
	$(CC) $(CFLAGS) $^ -o $@

# poisson-fem 的构建规则
poisson-fem: poisson-fem.o
	$(CC) $(CFLAGS) $^ -o $@

# learn 的构建规则（C++ 程序）
learn: learn.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# 通用规则：从 .c 文件生成 .o 文件
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# 通用规则：从 .cpp 文件生成 .o 文件
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理生成的文件
clean:
	rm -f *.o $(TARGETS)

# 伪目标声明
.PHONY: all clean
