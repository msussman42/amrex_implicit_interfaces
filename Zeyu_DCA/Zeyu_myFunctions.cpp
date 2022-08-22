#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;
#include "Zeyu_structfordrop.h"
#include "Zeyu_myFunctions.h"

void Initialization(vector<dropCondition> &droplet)
{
    string line;
    string paramName;
    double paramValue;

    ifstream input ("input_date_file_Zeyu");
    if(input.is_open()){
        while(getline(input, line)){
            if(line[0] != '#'){
                stringstream buffer;
                buffer << line;
                buffer >> paramName;
                if(paramName.compare("drop") == 0){
                    dropCondition temp;
                    buffer >> paramValue;
                    temp.x = paramValue;
                    buffer >> paramValue;
                    temp.y = paramValue;
                    buffer >> paramValue;
                    temp.z = paramValue;
                    buffer >> paramValue;
                    temp.r = paramValue;
                    droplet.push_back(temp); 
                }
                else
                    cerr << "unknown parameter: " << paramName << endl;
            }
        }
    input.close();
    }
    else
        cerr << "input file not found" << endl;
}

void Initial_Level_Set_Function_Thermal_Spray(const double x, const double y, const double z, double &phi, const vector<dropCondition> droplet)
{
    int in_droplet = 0;
    double min = 1e20;
    double temp = 0.0;
    for(unsigned int n = 0; n < droplet.size(); ++n){
        temp = droplet[n].r-sqrt(pow(x-droplet[n].x,2)+pow(y-droplet[n].y,2)+pow(z-droplet[n].z,2));
        if(temp > 0.0){
            phi = temp;
            in_droplet = 1;
            break;
        }
        else if(std::abs(temp) < min)
            min = std::abs(temp);
    }
    if(in_droplet == 0) phi = -min;
}
