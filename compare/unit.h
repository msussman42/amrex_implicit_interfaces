#ifndef  _UNIT_
#define  _UNIT_
#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>

class unit {
	public:
	double x_vel,y_vel;  //velocity for liquid
	double X_vel,Y_vel,phi_vapor;   //velocity for vapor


	public:
	unit( double a,double b,double c,double d,double e );
	unit()  ;
	unit (const unit &);

	void setValue( double a,double b,double c,double d,double e);
	~unit() {};
	
	friend unit operator + (const unit &p,const  unit &q) ;
	friend unit operator - (const unit &p,const  unit &q) ;
	friend unit operator * (const unit &p, const  unit &q) ;
	friend unit operator * (const unit p, double f) ;
	friend unit operator / (const unit &p, const unit &q) ;
//	unit & operator = (const unit &q) ; 
	friend ostream & operator << (ostream &, unit &);
};

#endif



	
