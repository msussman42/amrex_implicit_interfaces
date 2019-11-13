#include <iostream>
#include <vector>
using namespace std;

#include "Zeyu_struct.h"
#include "Zeyu_myFunctions.h"

/*

CC = g++
CPPFLAGS = -std=c++11
OPT = -O3 -Wall

SRC = Zeyu_main.cpp
OBJ = Zeyu_main.o

all: main
	@echo "----Compilation Done----"

main:
	$(CC) $(CPPFLAGS) $(OPT) Zeyu_main.cpp Zeyu_myFunctions.cpp -o Zeyu_main

clean:
	rm -f *.o Zeyu_main
 
 */

int main()//unit of length is [mu.m]
{
    #include "Zeyu_variables.h"

    Initialization(l, n, ndrop, droplet);//read input file

    Grid(xc, yc, zc, l, n);//set up mesh

    Initial_Level_Set_Function_Thermal_Spray(xc, yc, zc, phi, droplet);//calculate signed distance function

    Output(xc, yc, zc, phi);//output

    return 0;
}
