#ifndef   _DATASET_
#define   _DATASET_
#include "unit.h"   

class dataSet
{
	private:
		int M,N;
		unit *data;
	public:
		dataSet() {data = NULL; };
		dataSet(int,int);
		dataSet(char  *filename );
		dataSet(dataSet &v);

		~dataSet() {delete [] data;  };//have to use [] for class objects
		dataSet  refine( );
		dataSet  coarsen();
		double   getError( dataSet & q, int item);

		double interface_error(dataSet &q,int areaflag);
		double velocity_error(dataSet &q,int phase);
		double error(dataSet &q, int item);
		
                void Interface(dataSet &p, dataSet &q);



		int getM() {return M;}
		int getN() {return N;}

		friend ostream & operator << (ostream &, dataSet &);
		friend dataSet operator + ( const dataSet &, const dataSet &);
		friend dataSet operator - ( const dataSet &, const dataSet &);
		unit& operator ()  ( int ,int ) ;

		dataSet& operator = (const dataSet &);
};

#endif
