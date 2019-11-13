#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>
using namespace std;
#include "Zeyu_struct.h"
#include "Zeyu_myFunctions.h"

void Initialization(double *l, int *n, int &ndrop, vector<dropCondition> &droplet)
{
    string line;
    string paramName;
    double paramValueDouble;
    int    paramValueInt;

    ifstream input ("inputs_data_file_Zeyu");
    if(input.is_open()){
        while(getline(input, line)){
            if(line[0] != '#'){
                stringstream buffer;
                buffer << line;
                buffer >> paramName;
                if(paramName.compare("l") == 0){
                    buffer >> paramValueDouble;
                    l[0] = paramValueDouble;
                    buffer >> paramValueDouble;
                    l[1] = paramValueDouble;
                    buffer >> paramValueDouble;
                    l[2] = paramValueDouble;
                }
                else if(paramName.compare("n") == 0){
                    buffer >> paramValueInt;
                    n[0] = paramValueInt;
                    buffer >> paramValueInt;
                    n[1] = paramValueInt;
                    buffer >> paramValueInt;
                    n[2] = paramValueInt;
                }
                else if(paramName.compare("ndrop") == 0){
                    buffer >> paramValueInt;
                    ndrop = paramValueInt;
                }
                else if(paramName.compare("drop") == 0){
                    dropCondition temp;
                    buffer >> paramValueDouble;
                    temp.x = paramValueDouble;
                    buffer >> paramValueDouble;
                    temp.y = paramValueDouble;
                    buffer >> paramValueDouble;
                    temp.z = paramValueDouble;
                    buffer >> paramValueDouble;
                    temp.r = paramValueDouble;
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

void Grid(vector<double> &xc, vector<double> &yc, vector<double> &zc, const double *l, const int *n)
{
    // x-coordinates
    xc.push_back(0.0);
    for(int i = 1; i <= n[0]; ++i)
        xc.push_back((i-0.5)*l[0]/n[0]);
    xc.push_back(l[0]);
    // y-coordinates
    yc.push_back(0.0);
    for(int i = 1; i <= n[1]; ++i)
        yc.push_back((i-0.5)*l[1]/n[1]);
    yc.push_back(l[1]);
    // z-coordinates
    zc.push_back(0.0);
    for(int i = 1; i <= n[2]; ++i)
        zc.push_back((i-0.5)*l[2]/n[2]);
    zc.push_back(l[2]);
}

void Initial_Level_Set_Function_Thermal_Spray(const vector<double> xc, const vector<double> yc, const vector<double> zc, vector<double> &phi,const vector<dropCondition> droplet)
{
    vector<vector<double>> phiforall;//calucate phi for each droplet respectively.
    for(unsigned int n = 0; n < droplet.size(); ++n){
        vector<double> phiforeach;
        for(unsigned int i = 0; i < xc.size(); ++i){
            for(unsigned int j = 0; j < yc.size(); ++j){
                for(unsigned int k = 0; k < zc.size(); ++k){
                    phiforeach.push_back(droplet[n].r-sqrt(pow(xc[i]-droplet[n].x,2)+pow(yc[j]-droplet[n].y,2)+pow(zc[k]-droplet[n].z,2)));
                }
            }
        }
        phiforall.push_back(phiforeach);
    }
    //calucate phi for whole field
    for(unsigned int i = 0; i < xc.size(); ++i){
        for(unsigned int j = 0; j < yc.size(); ++j){
            for(unsigned int k = 0; k < zc.size(); ++k){
                int in_droplet = 0;
                double min = 1e20;
                for(unsigned int n = 0; n < droplet.size(); ++n){
                    if(phiforall[n][i*yc.size()*zc.size()+j*zc.size()+k] > 0.0){
                        phi.push_back(phiforall[n][i*(yc.size()*zc.size())+j*zc.size()+k]);
                        in_droplet = 1;
                        break;
                    }
                    else if(abs(phiforall[n][i*yc.size()*zc.size()+j*zc.size()+k]) < min)
                        min = abs(phiforall[n][i*yc.size()*zc.size()+j*zc.size()+k]);
                }
                if(in_droplet == 0) phi.push_back(-min);
            }
        }
    }
}

void Output(const vector<double> xc, const vector<double> yc, const vector<double> zc, const vector<double> phi)
{
    ofstream plot("plot.plt");
    
    if(plot.is_open()){
        plot << "title = \"initialize phi\"\n";
        plot << "variables = \"x\", \"y\", \"z\", \"phi\"\n";
        plot << "zone t = \"00000\"\ni = " << xc.size() << ", j = "<< yc.size() << ", k = " << zc.size() << ", f = point\n";
        for(unsigned int k = 0; k < zc.size(); ++k){
            for(unsigned int j = 0; j < yc.size(); ++j){
                for(unsigned int i = 0; i < xc.size(); ++i){
                    plot << scientific << setprecision(5) << xc[i] << ' ' << yc[j] << ' ' 
                         << zc[k] << ' ' << phi[i*yc.size()*zc.size()+j*zc.size()+k] << '\n';
                }
            }
        }
        plot.close();
    }
    else
        cerr << "input file not found" << endl;
}
