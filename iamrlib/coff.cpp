//
//
// COnvert Fortran77 files to Fortran95
//
// Alan Kuhnle
// 5-5-10

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

void splits( string& s_currl, string& s_nexl, size_t spos ) {
  s_nexl = s_currl.substr( spos, string::npos );
  // Take care of indentation
  s_nexl.insert( 0, s_currl , 0, s_currl.find_first_not_of( " \t" ) );
  s_currl = s_currl.substr( 0, spos );
  //Add continuation character to current line
  s_currl += " &";  

  return;
}

int main(int argc, char** argv) 
{
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " input file" << endl;
    exit(1);
  }

  string ifname( argv[1] );
  string ofname = ifname + "95";

  ifstream if1( ifname.c_str() );
  ofstream of1( ofname.c_str() );

  string sline;
  bool b_cont;
  string s_currl, s_nexl;
  const string s_rep0( "probdata.H" );
  const string s_rep1( "probdataf95.H" );
  const string s_mod0( "DIMS" );

  b_cont = getline( if1, s_nexl );
  while ( b_cont ) {
    // Current line = next line
    s_currl = s_nexl;

    if (s_currl[0] == 'c') 
      s_currl[0] = '!';

    size_t spos = s_currl.find( s_rep0 );
    if (spos != string::npos) {
      s_currl.replace(spos, s_rep0.length(), s_rep1);
    } 

    spos = s_currl.find( s_mod0 );
    while (spos != string::npos) {
      // Need to place DIMS on its own line
      splits( s_currl, s_nexl, spos );

      //Write current line
      of1 << s_currl << endl;

      s_currl = s_nexl;

      //search for end of DIMS macro
      spos = s_currl.find( ')' );
      if (spos == string::npos) {
	//this shouldn't happen
	cerr << "Found end of string when searching for ')'\n";
	cerr << "Error.\n";
	exit( 2 );
      }
      
      if (s_currl.length() > (spos + 1)) {
	splits( s_currl, s_nexl, spos + 1 );

	of1 << s_currl << endl;
	s_currl = s_nexl;
      }
      else {
	//DIMS macro was the last thing on this line
	//No further division necessary

      }

      //DIMS macro is now on its own line
      //Search for any other DIMS macros
      spos = s_currl.find( s_mod0 );
    }
      

    b_cont = getline( if1, s_nexl );
    if ( b_cont ) {
      if (s_nexl.length() > 5 && s_nexl[5] != ' ' && s_nexl[0] == ' ') {
	s_currl += " &";
	s_nexl[5] = ' ';
      }
    }

    //Write current line to output file
    of1 << s_currl << endl;
  }

  if1.close();
  of1.close();

  return 0;
}



    



  
  
