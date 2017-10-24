// classes example
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <unistd.h>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <math.h>       /* log */
#include <ctime>
#include <algorithm> // use abs() needed for sort
#include <random>

//#include <string> // 
//#include<math.h> // use sqrt()



#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


using namespace std;


///
/// This program calculates the distribution of: 
///
/// (1) Number of first neighbours
/// (2) First neighbours distances
///
/// Based on the connectivity criterion imposed by topological clustering and NOT
/// on the distance based clustering
///

// Compilation g++ -std=c++11 Add_more_water.cpp	




double Uhc (int a)
{
	double Sigma = 2.6 ;
	unsigned int Grid = 500 ;
	double Lr= 20.0 ; //Angstroms
	double dr= Lr/Grid ; //Angstroms
	double r ;

	if ( a*dr >  Sigma ) { r = 0.0 ; }
	if ( a*dr <= Sigma ) { r = 1.0/Sigma ; }
	return r ;
}


double Gauss (int a)
{
	double Sigma = 2.6 ;
	unsigned int Grid = 500 ;
	double Lr= 20.0 ; //Angstroms
	double dr= Lr/Grid ; //Angstroms
	double r ;
	double Pi = 3.14159265359 ;
	r = 2.0/(Sigma*sqrt(2.0*Pi))*exp(-0.5*pow(a*dr/Sigma,2)) ;

	return r ;
}

//**********************************************************************************
// main function
//**********************************************************************************
int main( int argc, const char* argv[] )
{
	std::cout << endl ;
	std::cout << BOLDYELLOW << "    _________________________________________________" << std::endl ;
	std::cout << BOLDYELLOW << "            _=_                                      " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          q(-_-)p                                    " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          '_) (_`         Classical DFT for LJ       " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          /__/  \\         Carles Triguero 2017      " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "        _(<_   / )_                                  " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "       (__\\_\\_|_/__)                               " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "    _________________________________________________" << RESET << std::endl ;
	std::cout << BOLDYELLOW << "                                                     " << RESET << std::endl ;


	std::cout << endl ;

	// Parameters:
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------	
	#define kB 1
	#define Lambda  1
	#define Temperature 1
	#define h 0.0001
	#define Mu -10.0  //-3.41125
	#define Epsilon 1.6		// kJ/mol
	#define Sigma 2.6		// Angstrom

	unsigned int Grid = 500 ;		// Number of atoms 60A/rm
	unsigned int Loop = 10000 ;	// Minimization loop
	
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------		
	
	// Core
	double Lr= 20.0 ; //Angstroms
	double dr= Lr/Grid ; //Angstroms

	ofstream CoreOut ;

	CoreOut.open ("GaussCore.dat", ios::out | ios::trunc);
	for ( int k = 0; k < Grid; ++k ) { CoreOut << k*dr << "\t" << Gauss(k) << std::endl ; }
	CoreOut.close ();
	
	CoreOut.open ("HardCore.dat", ios::out | ios::trunc);
	for ( int k = 0; k < Grid; ++k ) { CoreOut << k*dr << "\t" << Uhc(k) << std::endl ; }
	CoreOut.close ();



	
	for ( int k = 1; k < argc ; ++k )
	{
		if ( ( argv[k] == std::string("-h") ) || ( argv[k] == std::string("-HELP") ) || ( argv[k] == std::string("-H") ) || ( argv[k] == std::string("-help") ) )
		{ 
			cout  << BOLDBLACK << "    HELP:" << RESET << std::endl ;
			cout << "    Generates a xyz configuration of A and B species" << std::endl ;
			cout << "    Interpenetrated rectangular 2d lattices" << std::endl ;
			cout << "    Lattice parameter is 1" << std::endl ;
			cout << "    Lattice parameter is 1" << std::endl ;
			cout << "    X-+-X-+-" << std::endl ;
			cout << "    +-O-+-O-" << std::endl ;
			cout << "    X-+-X-+-" << std::endl ;
			cout << "    ------------------------------------------------" << std::endl ;
			cout << BOLDBLACK << "    Mandatory flags:" << RESET << std::endl ; 
			cout << "    -m" << "\t" << "    to set molecule number (e.g. ./a.out -m 1000)" << std::endl ;
			cout << BOLDBLACK << "    Optional flags:" << RESET << std::endl ;
			cout << "    -h" << "\t" << "    to get help  (e.g. ./a.out -help)" << std::endl ; 
			cout << std::endl ;
			return (0) ;
		}
		//if ( ( argv[k] == std::string("-lx") ) || ( argv[k] == std::string("-Lx") )  || ( argv[k] == std::string("-LX") ) ) { Lx = atoi(argv[k+1]) ; }
		//if ( ( argv[k] == std::string("-ly") ) || ( argv[k] == std::string("-Ly") )  || ( argv[k] == std::string("-LY") ) ) { Ly = atoi(argv[k+1]) ; }
	}



	// Count execution time
	int start_s=clock();
	
	// Files to store data
	ofstream DensityOut ;
	DensityOut.open ("Density.dat", ios::out | ios::trunc);

	// Initial density
	std::vector<double> Rho ;
	for ( int k = 0; k < Grid; ++k ) { Rho.push_back (1.0) ; }

	double U ;
	double beta = 1/(kB*Temperature) ;
	double K1, K2, K3, K4 ;

	unsigned int count = 0 ;
	for ( int i = 0; i < Loop; ++i ) // Minimization Loop
	{
		for ( int k = 0; k < Grid; ++k ) // Spatial Loop
		{
			double y = Rho[k] ;

			//U = Uhc (k) ;
			U = Gauss (k) ;

//			t=k
//			y=Rho

			K1 = 1.0/beta*log(pow(Lambda,3)*y)+U-Mu ;
			K2 = 1.0/beta*log(pow(Lambda,3)*(y+h/2.0*K1))+U-Mu ;
			K3 = 1.0/beta*log(pow(Lambda,3)*(y+h/2.0*K2))+U-Mu ;
			K4 = 1.0/beta*log(pow(Lambda,3)*(y+h*K3))+U-Mu ;
			Rho[k] = y + h/6.0*( K1 + 2.0*K2 + 2.0*K3 + K4 ) ;
		}
	}

	
	for ( int k = 0; k < Grid; ++k ) { DensityOut << k*dr << "\t" << Rho[k] << std::endl ; }

	


	// End counting time 	
 	int stop_s=clock();
 	
	// Execution time
 	std::cout << endl ;
	std::cout  << YELLOW << "    -> Execution time: " << BOLDYELLOW << ( stop_s - start_s )/double(CLOCKS_PER_SEC) << RESET << YELLOW << " seconds" << std::endl ;
	std::cout << endl ;
	std::cout  << BOLDYELLOW << "    Program finished" << RESET << std::endl ;
	std::cout << endl ;
	
	return (0) ;
}
