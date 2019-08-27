	//Program that reads in a data file, runs the data through a function, then outputs
	// the new data to a new file.

	#include <iostream>
	#include <cmath>
	#include <fstream>

	using namespace std;
	
	double data[2][2];
	double RHS[2][1];
	double A_inverse[2][2];
	double y_i[1][2];

	class Cramer
	{
		public:
			//Setting accessors and matrices
			void find_det();
			void inverse_mat();
			void get_adjugate();
			void get_SI();
		private:
			int i,j;
			double det;
			double SI_mat[2][1];
			double A_adjugate[2][2];
	};


	//Find the determinant of data[2][2]
	void Cramer::find_det()
	{
		//Go back to main program and create Linked list/nodes/pointers
		//to store the valus of data[][] at each point?
		
		det = (data[0][0]*data[1][1])-(data[0][1]*data[1][0]);
		return;

	}//end of get_det


	//Create adjugate matrix of data[2][2].
	void Cramer::get_adjugate()
	{
		//Initialize matrix.
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				A_adjugate[i][j]=0.0;
			}//end of for
		}//end of for

		//Assign the values of the adjugate matrix
		
		A_adjugate[0][0] = data[1][1];
		A_adjugate[0][1] = -1*data[0][1];
		A_adjugate[1][0] = -1*data[1][0];
		A_adjugate[1][1] = data[0][0];
		
		return;

	}//end of adjugate


	//Create inverse function
	void Cramer::inverse_mat()
	{
		//Intialize inverse matrix.
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				A_inverse[i][j]=0.0;
			}//end of i for
		}//end of j for

	
		get_adjugate();
		find_det();

		cout << "\nDet: \t" << det << endl;
		cout << "\n\nA_inverse is:\n";

		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				A_inverse[i][j] = (1/det)*A_adjugate[i][j];
			
				cout << A_inverse[i][j] << "\t";
			}//end of j for
			
			cout << endl;
		}//end of i for
		return;

	}//End of inverse_mat
	

	void Cramer::get_SI()
	{
		//Initialize SI_mat
		for(i=0;i<2;i++)
		{
			for(j=0;j<1;j++)
			{
				SI_mat[i][j] = 0.0;
			}//end of j for
		}//end of i for
		
		//Call inverse function.
		inverse_mat();

		cout << "\n\nThe Slope-Intercept matrix is: \n";

		//Calculate Slope-Intercept values & return.
		for(i=0;i<2;i++)
		{
   	         SI_mat[i][0] =( (A_inverse[i][0]*RHS[0][0]) + (A_inverse[i][1]*RHS[1][0]) );

                 if (i==0) {
                  cout << "intercept (log(A)) " << SI_mat[i][0] << endl;
                  cout << "intercept A " << exp(SI_mat[i][0]) << endl;
                 } else {
                  cout << "slope " << SI_mat[i][0] << endl;
                 }

		}//end of i for
		
		cout << endl;

		//Caclculate Total Error; use equation 8.5
		double Total_Error1 = pow( y_i[0][1]-(exp(SI_mat[0][0])*pow(data[1][0],SI_mat[1][0])), 2);
		cout << "Total Error using equation 8.5 is: " << "\t" << Total_Error1 << "\n" << endl;
		double Total_Error2 = pow( y_i[0][1]-(exp(SI_mat[0][0])*exp(SI_mat[1][0]*data[1][0])), 2);
		cout << "Total Error using equation 8.4 is: " << "\t" << Total_Error2 << "\n" << endl;
		return;
	}//end of get_SI



	int main()	
	{
		
		//Method for reading in a file, performing a function on it, 
		//and then writing it out to a new file.

		//User enters input data file.
		char infileName[20];
		cout << "\nEnter the name of the Data Input file: ";
		cin.get (infileName,20);
		cout << "\n";

		//User creates output data file
		char filename[50];
		cout << "Enter the name of the Output data file: ";
		cin >> filename; //Name of the file the datapoints are going to
		cout << "\n";

		ofstream fout(filename); //Opens filename for writing.



		//Create streamobject.
		ifstream infile; 
		infile.open(infileName);

		int n=2;
		int m=20000;
		int rows_read=0;

		//Pre-allocate memory for the matrix.
		double **mat = new double*[m];

		for(int i=0; i<m; i++)
		{
 			mat[i] = new double[n];
		}//End of for.


		//Create actual matrix
  		while (infile.good())
		{
			//following reads in the time for the rows_read line
			infile >> mat[rows_read][0];

			//following reads in the amplitude for rows_read line
			infile >> mat[rows_read][1];

			rows_read++;
			if (rows_read >= 20000)
			{
				cout << "Error has occurred.\n";
			}//end of if

		}//End of while.


		//Initiate data. Create for loops for i=0 to 1, j=0 to 1, data[i][j]=0
		int i,j;
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				data[i][j]=0.0;
			}//End of for.
		}//End of for.

		//Declare a 2x1 Array "RHS" for right-hand side, and Initialize.
		for(i=0;i<2;i++)
		{
			for(j=0;j<1;j++)
			{
				RHS[i][j]=0.0;
			}//End of for.
		}//End of for.

		//Intialize Error Array
		for(i=0;i<1;i++)
		{
			for(j=0;j<2;j++)
			{
				y_i[i][j]=0.0;
			}//End of for.
		}//End of for.

//look into log base 10 instead of ln

		//Create matrix for 'loop=0 to (rows_read-1)'
		int loop;
		for(loop=0; loop < (rows_read-1); loop++)
		{
			//data [0][0] is incremented by 1--> 
			data[0][0] = data[0][0] + 1;
			//data [0][1] is incremented by mat[loop][0] -->
			data[0][1] = data[0][1] + mat[loop][0];
			//data [1][0] is incremented by mat[loop][0] --> 
			data[1][0] = data[1][0] + mat[loop][0];
			//data [1][1] is incremented by mat[loop][0]^2 --> 
			data [1][1] = data[1][1] + pow(mat[loop][0],2);
					
			//RHS [0][0] is incremented by log (mat[loop][1])
			RHS[0][0] = RHS[0][0] + log(mat[loop][1]);
			//RHS [1][0] is incremented by mat[loop][0]*log(mat[loop][1])
			RHS[1][0] = RHS[1][0] + mat[loop][0]*log(mat[loop][1]);	

			y_i[0][1] = y_i[0][1] + mat[loop][1];		
		}//End of 'loop' for.

		cout << "\nThe Data matrix is: " << "\n";

		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				cout << data[i][j] << "\t";
			}//end of j for
			
			cout << endl;
		}//end of i for.

		cout << "\n\nThe Right-Hand Side matrix is: " << "\n";

		for(i=0;i<2;i++)
		{
			for(j=0;j<1;j++)
			{
				cout << RHS[i][j] << "\t";
			}//end of j for
			
			cout << endl;
		}//end of i for.


		//Use cramer's formula to create inverse matrix. = (1/det)*A_t = A^-1	
		Cramer alpha;
		alpha.find_det();
		alpha.get_adjugate();

		//Multiply by RHS---this provides the slope and the intercept
		//for the best fit for the log data.
		alpha.get_SI();	


		double t_v = 1.0;

		//Declare new array N_v and initialize
		double N_v[1][1];
		for(i=0;i<1;i++)
		{
			for(j=0;j<1;j++)
			{
				N_v[i][j] = 0.0;
			}//end of j for
		}//end of i for.

		//Create loop that calculates the value of Nv located within
		//the linear regime of the graph "log(RHS[0][1])" vs "data[1][0]/t_v"
		//and write these points to an output file
		for(int loop2=0; loop2 < (rows_read-1); loop2++)
		{
			//N_v[0][0] is incremented by 1--> 
			N_v[0][0] = (mat[loop2][0]/t_v);
			fout << N_v[0][0] << "\t";
			//N_v[0][1] is incremented by log(mat[loop][0]) -->
			N_v[0][1] = log(mat[loop2][1]);
			fout << N_v[0][1] << endl;
		}//End of loop2 for loop.


		//create another instance of Cramer that:
			//Creates a loop that plots the points of "N_v" vs data[1][0]/t_v"
			//and write these points to an output file


		//Close output file
		fout.close();

		infile.close();
	
		return 0;

	}//End of main.
