
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace std;
ofstream ofline;


// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

// Function to read in data from screen
void read_input(int&, int&, double&, double&,double&);

// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);

// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);

// prints to file the results of the calculations
void output(int, int, double, double *);


// main program
int main()
{

    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1){
        cout <<"Bad Usage: " << argv[0]<< "read also output file on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }

    ofile.open(outfilename);

    //    Read in initial values such as size of lattice, temp and cycles
    read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1; // random starting point

    for ( double temp = initial_temp; temp <= final_temp; temp += temp_step){

        //initialize energy and magnetization
        E = M = 0;

        // setup array for possible energy changes
        for( int de = -8; de <= 8; de++ ) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        // initialize array for expectation values
        for (int i = 0; i < 5; i++) average[i] = 0;
        initialize(n_spins, double temp, spin_matrix, E, M);

        // start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w);

            // update expectation values
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);

        }
        // print results
        output(n_spins, mcs, temp, average);

    }
    free_matrix((void **) spin_matrix); // free memory
    ofile.close(); // close output file
    return 0;

}

