CC = g++
CFLAGS = -Wall -Wextra -std=c++11 -g -w -O2 -fopenmp

SRCS = radio.cpp Graph.cpp ./BuildIndex/TreeIndex.cpp ./BuildIndex/hierarchy.cpp ./Algorithm/Online.cpp ./BuildIndex/LCS.cpp ./core-maint/gadget/gadget.cc ./core-maint/gadget/treap.cc ./core-maint/gadget/heap.cc ./core-maint/glist/glist.cc
HDRS = Graph.h CoreGroup.h TreeIndex.h Online.h hierarchy.h LCS.h glist.h treap.h heap.h gadget.h

OBJS = $(SRCS:.cpp=.o) $(SRCS:.cc=.o)  # 正确处理 .cpp 和 .cc 文件
OBJS := $(filter %.o, $(OBJS))  # 确保只包含 .o 文件


TARGET = a.out

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(TARGET)
