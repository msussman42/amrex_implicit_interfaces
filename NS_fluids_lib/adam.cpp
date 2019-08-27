#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <iomanip>
#include <ctime>

// Alex Ayoub October, 2018

using namespace std;

// Linear Algebra Vector Operations
double *addscalar(double c, double x[], int n) // this looks like (constant+vector)
{
	double *z = new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = c+x[i];
	}
	return z;
}

double *addrootscalar(double c, double x[], int n) // this looks like (vector + constant)^1/2
{
	double *z = new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = sqrt(c+x[i]);
	}
	return z;
}

double *scalar(double c, double x[], int n)// this looks like (constant*vector)
{
	double *z = new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = c*x[i];
	}
	return z;
}

double *addvec(double x[], double y[], int n)// this looks like (vector+vector)
{
	double *z=new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

double *elementwisesquare(double x[], double y[], int n)// this looks like (vector[i]*vector[i])
{
	double *z=new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = x[i] * y[i];
	}
	return z;
}

double *divideVec(double x[], double y[], int n)// this looks like (vector[i]/vector[i])
{
	double *z=new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = x[i] / y[i];
	}
	return z;
}


//ADAM Begins here:



double *f(double x[], int n)
{
	double *z = new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] =  x[i]*x[i] - 4.0*x[i] + 4.0; //function we are trying to optimize, should converge to x=2
	}
	return z;
}

double *fprime(double x[], int n)
{
	double *z = new double[n];
	for(int i=0;i<n;i++)
	{
		z[i] = 2.0*x[i] - 4.0;  //where to add the gradient of our lost(cost function), also should make a seperate program(.cpp file) to do complex step differentiation or automatic differentiation 
	}
	return z;
}


// This routine finds the minimum of f(Theta) using fprime(Theta).
// inputs:
//  Theta[] is an initial guess
//  alpha is the step size
//  beta_1, beta_2 are exponential decays rates:
//    e.g. beta_1=0.9 and beta_2=0.999  alpha=0.01
double *adam(double Theta[], double alpha, double beta_1, double beta_2, double eps, int n)
{

	//initializing the 1st and 2nd moment vectors to be 0.0
	double * m_t = new double[n];
	double * v_t = new double[n];
	double * theta = new double[n];
	double * theta_prev = new double[n];
	double * g_t = new double[n];
	double * m_tHat = new double[n];
	double * v_tHat = new double[n];
	for(int i=0;i<n;i++)
	{
		m_t[i] = 0.0;
		v_t[i] = 0.0;
		theta[i] = Theta[i];
	}
	
	//initialize time step
	int t = 0;
	
	//initialize convergence check
	int convergeCheck = 0;
	
	while(t<7000) //need to decide what our critera for convergence should be, right now the program just iterates through 100000 time steps
	{
		t =t + 1; 
		g_t = fprime(theta,n); //gets the gradient of the "stochastic function", we need to either use complex step derivative or reserve mode automatic differientation or just compute fprime
		m_t = addvec(scalar(beta_1,m_t,n),scalar(1.0-beta_1,g_t,n),n); //update moving averages of the gradients
		v_t = addvec(scalar(beta_2,v_t,n),scalar(1.0-beta_2,elementwisesquare(g_t,g_t,n),n),n); //updates the moving average of the gradient squared
		m_tHat = scalar(pow(1.0-pow(beta_1,t),-1.0),m_t,n); // calculates the bias-correct estimates
		v_tHat = scalar(pow(1.0-pow(beta_2,t),-1.0),v_t,n); // calculates the bias-correct estimates
		theta_prev = theta;
		theta = addvec(theta,divideVec(scalar(-1.0*alpha,m_tHat,n),addrootscalar(eps,v_tHat,n),n),n); //updates the parameters
	}
	return theta; //returns the minimum of the function
}

int main()
{
	// initialize parameters to input into our adam function
	const int n = 10;
	double Theta[n] = {1.0,3.0,3.4444,25.66,-43.55,6.0,-12.3,.0098,0.0,-10.0};
	double * theta = new double[n];
	theta = Theta;
	double alpha = 0.01;
	double beta_1 = 0.9;
	double beta_2 = 0.999;
	double eps = pow(10.0,-8.0);
	int start_s=clock();
	theta = adam(theta,alpha,beta_1,beta_2,eps,n);
	int stop_s=clock();
	cout << "time(milliseconds): " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;

	for(int i=0;i<n;i++)
	{
		cout << abs(2.0 - theta[i]) << endl; //prints the error of the adam optimizer output.
	}
}

