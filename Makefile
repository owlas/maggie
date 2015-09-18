#
# Make file for compiling the maggie executable
#

CC=g++
INC_PATH=include
OBJ_PATH=objects
LIB_PATH=lib
SRC_PATH=src

GTEST_DIR=gtest-1.7.0

LIB_DIR = .
LIBS=-lgtest -lpthread -lmaggie

CPP_FLAGS=--std=c++11 -W -Wall -pedantic

SOURCES=$(wildcard $(LIB_PATH)/*.cpp)
OBJ_FILES=$(addprefix $(OBJ_PATH)/,$(notdir $(SOURCES:.cpp=.o)))

# Default invokes the primary target
#
default: main

# gtest builds the gtest module
#
gtest:
	$(CC) -I$(GTEST_DIR)/include -I$(GTEST_DIR) -pthread -c $(GTEST_DIR)/src/gtest-all.cc
	ar -rv libgtest.a gtest-all.o

# Create the main executable - main
#
main: src/main.cpp libmaggie.so
	$(CC) $(CPP_FLAGS) -I$(INC_PATH) -I$(GTEST_DIR)/include -o $@ $(SRC_PATH)/main.cpp -L$(LIB_DIR) $(LIBS)

# Shared library maggie.so used for objects
#
libmaggie.so: $(OBJ_FILES)
	$(CC) -shared -o $@ $^



# Compile objects from source
#
$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CC) $(CPP_FLAGS) -I$(INC_PATH) -c -fPIC -o $@ $<

# Clean up executable files
#
clean:
	rm -f objects/*
	rm -f main

cleangtest:
	rm -f gtest-all.o, libgtest.a
