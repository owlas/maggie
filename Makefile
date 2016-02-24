#
# Make file for compiling the maggie executable
#

CC=g++
INC_PATH=include
OBJ_PATH=objects
LIB_PATH=lib
SRC_PATH=src

GTEST_DIR=googletest/googletest
GMOCK_DIR=googletest/googlemock

LIB_DIR = .
LIBS=-lgmock -lpthread -lmaggie

CPP_FLAGS=--std=c++11 -W -Wall -pedantic -Wno-unused-local-typedefs -g

SOURCES=$(wildcard $(LIB_PATH)/*.cpp)
OBJ_FILES=$(addprefix $(OBJ_PATH)/,$(notdir $(SOURCES:.cpp=.o)))

# Default invokes the primary target
#
default: runtests convergence_tests llg_solver_convergence llgsim


initial:
# Create objects folder and pull inih submodule
	mkdir -p objects
	git submodule init
	git submodule update
	
# gtest builds the gtest and gmock modules
#
gmock:
	$(CC) 	-isystem $(GTEST_DIR) -I$(GTEST_DIR)/include \
		-isystem $(GMOCK_DIR) -I$(GMOCK_DIR)/include \
		-pthread -c $(GTEST_DIR)/src/gtest-all.cc
	$(CC) 	-isystem $(GTEST_DIR) -I$(GTEST_DIR)/include \
		-isystem $(GMOCK_DIR) -I$(GMOCK_DIR)/include \
		-pthread -c $(GMOCK_DIR)/src/gmock-all.cc
	ar -rv libgmock.a gtest-all.o gmock-all.o

# Create the tests executable - tests
#
runtests: src/tests.cpp libmaggie.so
	$(CC) 	$(CPP_FLAGS) \
		-I$(INC_PATH) -I$(GTEST_DIR)/include -I$(GMOCK_DIR)/include \
		-o $@ $(SRC_PATH)/tests.cpp -L$(LIB_DIR) $(LIBS)

# Create the convergence tests
#
convergence_tests: src/convergence_tests.cpp libmaggie.so
	$(CC) 	$(CPP_FLAGS) \
		-I$(INC_PATH) -I$(GTEST_DIR)/include -I$(GMOCK_DIR)/include \
		-o $@ $(SRC_PATH)/convergence_tests.cpp -L$(LIB_DIR) $(LIBS)

# Create the llg convergence tests
#
llg_solver_convergence: src/llg_solver_convergence.cpp libmaggie.so
	$(CC) 	$(CPP_FLAGS) \
		-I$(INC_PATH) -I$(GTEST_DIR)/include -I$(GMOCK_DIR)/include \
		-o $@ $(SRC_PATH)/llg_solver_convergence.cpp -L$(LIB_DIR) $(LIBS)


# Create the examples
#
llgsim: src/llgsim.cpp libmaggie.so
	$(CC) 	$(CPP_FLAGS) \
		-I$(INC_PATH) -I$(GTEST_DIR)/include -I$(GMOCK_DIR)/include \
		-o $@ $(SRC_PATH)/llgsim.cpp -L$(LIB_DIR) $(LIBS)

# Shared library maggie.so used for objects
#
libmaggie.so: $(OBJ_FILES)
	$(CC) -shared -o $@ $^



# Compile objects from source
#
$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CC) 	$(CPP_FLAGS) \
		-I$(INC_PATH) -c -fPIC \
		-o $@ $<

# Clean up executable files
#
clean:
	rm -f objects/*
	rm -f main

gclean:
	rm -f gtest-all.o, libgtest.a
