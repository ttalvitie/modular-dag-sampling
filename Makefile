CXX ?= g++
CFLAGS ?= -O2 -Wall -Wextra -pedantic -march=native -std=c++11
LDFLAGS ?= 
COMMONSRCS := $(shell find src -name '*.cpp' -not -path 'src/sampler.cpp')
COMMONOBJS := $(COMMONSRCS:%.cpp=%.o)
SRCS = $(COMMONSRCS) src/sampler.cpp
OBJS := $(SRCS:%.cpp=%.o)
DEPS := $(SRCS:%.cpp=%.d)

.PHONY: all

all: sampler

sampler: $(COMMONOBJS) src/sampler.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -f sampler $(OBJS) $(DEPS)

-include $(DEPS)
