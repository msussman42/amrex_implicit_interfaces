#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
#include "Zeyu_read_input_file_dynamic_contact_angle.h"

void read_input_file(int &dim, double &mu_l, double &mu_g, double &sigma, double &thet_s, int &imodel, int &ifgnbc, int &idelta, int &ifmicro, vector<double> &param)
{
    string line;
    string paramName;
    double paramReal;
    int paramInt;

    //defualt value
    imodel = 1;
    ifgnbc = 0;
    idelta = 1;
    ifmicro = 1;

    ifstream input ("input_date_file_Zeyu_dynamic_contact_angle");
    if(input.is_open()){
        while(getline(input, line)){
            if(line[0] != '#' && !line.empty()){
                stringstream buffer;
                buffer << line;
                buffer >> paramName;
                if(paramName.compare("dim") == 0){
                    buffer >> paramInt;
                    dim = paramInt;
                }
                else if(paramName.compare("mu_l") == 0){
                    buffer >> paramReal;
                    mu_l = paramReal;
                }
                else if(paramName.compare("mu_g") == 0){
                    buffer >> paramReal;
                    mu_g = paramReal;
                }
                else if(paramName.compare("sigma") == 0){
                    buffer >> paramReal;
                    sigma = paramReal;
                }
                else if(paramName.compare("thet_s") == 0){
                    buffer >> paramReal;
                    thet_s = paramReal;
                }
                else if(paramName.compare("imodel") == 0){
                    buffer >> paramInt;
                    imodel = paramInt;
                }
                else if(paramName.compare("ifgnbc") == 0){
                    buffer >> paramInt;
                    if(imodel != 1)
                        ifgnbc = 0;
                    else
                        ifgnbc = paramInt;
                }
                else if(paramName.compare("idelta") == 0 && ifgnbc != 0){
                    buffer >> paramInt;
                    idelta = paramInt;
                }
                else if(paramName.compare("ifmicro") == 0 && ifgnbc != 0){
                    buffer >> paramInt;
                    idelta = paramInt;
                }
                else if(paramName.compare("lamda") == 0 && imodel == 1){
                    buffer >> paramReal;
                    param.push_back(paramReal);
                }
                else if(paramName.compare("l_macro") == 0 && imodel == 1){
                    buffer >> paramReal;
                    param.push_back(paramReal);
                }
                else if(paramName.compare("l_micro") == 0 && imodel == 1){
                    buffer >> paramReal;
                    param.push_back(paramReal);
                }
                else if(paramName.compare("sigma_0") == 0 && imodel == 7){
                    buffer >> paramReal;
                    param.push_back(paramReal);
                }
                else if(paramName.compare("v_0") == 0 && imodel == 7){
                    buffer >> paramReal;
                    param.push_back(paramReal);
                }
                else
                    cout << "unused parameter: " << paramName << endl;
            }
        }
        input.close();
    }
    else
        cerr << "input file not found" << endl;
}

