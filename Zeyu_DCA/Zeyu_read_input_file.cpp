#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
#include "Zeyu_read_input_file.h"

void read_input_file(int &dim, double &mu_l, double &sigma, double &thet_s, double &thet_d_grid, double &l_grid, double &l_macro, double &l_micro, double &chi, double &x_cl, double &y_cl, double &z_cl)
{
    string line;
    string paramName;
    double paramValue;
    int paramInt;

    ifstream input ("input_date_file_Zeyu");
    if(input.is_open()){
        while(getline(input, line)){
            if(line[0] != '#'){
                stringstream buffer;
                buffer << line;
                buffer >> paramName;
                if(paramName.compare("dim") == 0){
                    buffer >> paramInt;
                    dim = paramInt;
                }
                else if(paramName.compare("mu_l") == 0){
                    buffer >> paramValue;
                    mu_l = paramValue;
                }
                else if(paramName.compare("sigma") == 0){
                    buffer >> paramValue;
                    sigma = paramValue;
                }
                else if(paramName.compare("thet_s") == 0){
                    buffer >> paramValue;
                    thet_s = paramValue;
                }
                else if(paramName.compare("thet_d_grid") == 0){
                    buffer >> paramValue;
                    thet_d_grid = paramValue;
                }
                else if(paramName.compare("l_grid") == 0){
                    buffer >> paramValue;
                    l_grid = paramValue;
                }
                else if(paramName.compare("l_macro") == 0){
                    buffer >> paramValue;
                    l_macro = paramValue;
                }
                else if(paramName.compare("l_micro") == 0){
                    buffer >> paramValue;
                    l_micro = paramValue;
                }
                else if(paramName.compare("chi") == 0){
                    buffer >> paramValue;
                    chi = paramValue;
                }
                else if(paramName.compare("x_cl") == 0){
                    buffer >> paramValue;
                    x_cl = paramValue;
                }
                else if(paramName.compare("y_cl") == 0){
                    buffer >> paramValue;
                    y_cl = paramValue;
                }
                else if(paramName.compare("z_cl") == 0){
                    buffer >> paramValue;
                    z_cl = paramValue;
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

