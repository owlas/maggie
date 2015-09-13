#
# Make file for compiling the maggie executable
#

INC_PATH=include/
OBJ_PATH=objects/

CPP_FLAGS=--std=c++11 -W -Wall -pedantic

# Default invokes the primary target
#
default: main


# Create the main executable - main
#
main: src/main.cpp objects/LangevinEquation.o objects/StocLLG.o
	g++ $(CPP_FLAGS) -I$(INC_PATH) src/main.cpp $(OBJ_PATH)/LangevinEquation.o $(OBJ_PATH)/StocLLG.o -o main


# LangevinEquation.cpp class is compiled to object
#
objects/LangevinEquation.o: src/LangevinEquation.cpp include/LangevinEquation.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/LangevinEquation.cpp -o objects/LangevinEquation.o

# StocLLG.cpp class is compiled to object
#
objects/StocLLG.o: src/StocLLG.cpp include/StocLLG.hpp
	g++ $(CPP_FLAGS) -I$(INC_PATH) -c src/StocLLG.cpp -o objects/StocLLG.o

# Clean up executable files
#
clean:
	rm objects/*
