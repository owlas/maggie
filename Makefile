#
# Make file for compiling the maggie executable
#

INC_PATH=include
OBJ_PATH=objects

GTEST_DIR=gtest-1.7.0

LIB_DIR = .
LIBS=-lgtest -lpthread

CPP_FLAGS=--std=c++11 -W -Wall -pedantic

# Default invokes the primary target
#
default: main

# gtest builds the gtest module
#
gtest:
	g++ -I$(GTEST_DIR)/include -I$(GTEST_DIR) -pthread -c $(GTEST_DIR)/src/gtest-all.cc
	ar -rv libgtest.a gtest-all.o


# Create the main executable - main
#
main: src/main.cpp objects/LangevinEquation.o objects/StocLLG.o objects/Integrator.o objects/RK4.o
	g++ $(CPP_FLAGS) -I$(INC_PATH) -I$(GTEST_DIR)/include -o main src/main.cpp -L$(LIB_DIR) $(LIBS) $(OBJ_PATH)/LangevinEquation.o $(OBJ_PATH)/StocLLG.o $(OBJ_PATH)/Integrator.o $(OBJ_PATH)/RK4.o


# LangevinEquation.cpp class is compiled to object
#
objects/LangevinEquation.o: src/LangevinEquation.cpp include/LangevinEquation.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/LangevinEquation.cpp -o objects/LangevinEquation.o

# StocLLG.cpp class is compiled to object
#
objects/StocLLG.o: src/StocLLG.cpp include/StocLLG.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/StocLLG.cpp -o objects/StocLLG.o

# Integrator.cpp class is compiled to object
#
objects/Integrator.o: src/Integrator.cpp include/Integrator.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/Integrator.cpp -o objects/Integrator.o

# RK4.cpp class is compiled to object
#
objects/RK4.o: src/RK4.cpp include/RK4.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/RK4.cpp -o objects/RK4.o

# Clean up executable files
#
clean:
	rm objects/*
