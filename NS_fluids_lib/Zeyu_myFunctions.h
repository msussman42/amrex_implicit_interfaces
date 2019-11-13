#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

void Initialization(double *, int *, int &, vector<dropCondition> &);

void Grid(vector<double> &, vector<double> &, vector<double> &, const double *, const int *);

void Initial_Level_Set_Function_Thermal_Spray(const vector<double>, const vector<double>, const vector<double>, vector<double> &, const vector<dropCondition>);

void Output(const vector<double>, const vector<double>, const vector<double>, const vector<double>);

#endif
