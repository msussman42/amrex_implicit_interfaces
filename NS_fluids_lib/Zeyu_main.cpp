#include <iostream>
#include <vector>
using namespace std;

#include "Zeyu_structfordrop.h"
#include "Zeyu_myFunctions.h"

/*
    g++ -std=c++11 -O3 -Wall Zeyu_main.cpp Zeyu_myFunctions.cpp -o Zeyu_main
*/

int main()//unit of length is [mu.m]
{
    vector<dropCondition> droplet;
    double phi;

    //read input file
    Initialization(droplet);

    //input  x, y and z from keyboard, calculate signed distance and output phi.
    double x, y, z;
    cout << "Please input coordinates of a point (x, y, z):" << endl;
    cin >> x >> y >> z;
    Initial_Level_Set_Function_Thermal_Spray(x, y, z, phi, droplet);
    cout << "\nThe value of phi of this point is: " << phi << endl;

    return 0;
}
