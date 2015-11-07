/*
   Program to solve the two-dimensional Ising model
   with zero external field.
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include "mpi.h"

using namespace  std;

ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}
// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
void output(int, int, double, double *);

int main(int argc, char* argv[])
{
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs, my_rank, numprocs;
    double w[17], average[5], total_average[5], initial_temp, final_temp, E, M, temp_step;
    //double accepted = 0;
    //double acceptedmoves = 0;

    // Read in output file, abort if there are too few command-line arguments
    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0 && argc <= 1) {
      cout << "Bad Usage: " << argv[0] <<
        " read output file" << endl;
      exit(1);
    }
    if (my_rank == 0 && argc > 1) {
      outfilename=argv[1];
      ofile.open(outfilename);
    }

    n_spins = 160; mcs = 1000000;  initial_temp = 2.1; final_temp = 2.7; temp_step =0.01;
    /*
    Determine number of intervall which are used by all processes
    myloop_begin gives the starting point on process my_rank
    myloop_end gives the end point for summation on process my_rank
    */
    int no_intervalls = mcs/numprocs;
    int my_loop_begin = my_rank*no_intervalls + 1;
    int my_loop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( my_loop_end < mcs) ) my_loop_end = mcs;

    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //Allocate memory for spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));

    //int startcount = 300000;

    //double energies[mcs-startcount-1];
    idum = -1-my_rank; // random starting point
    for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
        cout << "Running calculation for temperature = " << temperature << endl;
        //    initialise energy and magnetization
        E = M = 0.;
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
        // initialise array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        for( int i = 0; i < 5; i++) total_average[i] = 0.;
        initialize(n_spins, temperature, spin_matrix, E, M);
        // start Monte Carlo computation
        //double counter = 0;
        for (int cycles = my_loop_begin; cycles <= my_loop_end; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w);
            // update expectation values
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
            //counter += 1;
            //acceptedmoves += accepted;
            //if(cycles > startcount){
                //energies[cycles-startcount-1] = E;
                //ofile << setw(15) << energies[cycles-startcount-1] << endl;
            //}
            //if (counter == 1000){
            //ofile << setw(15) << setprecision(8) << cycles;
            //output(n_spins, cycles, temperature, average);
            //counter = 0;
            //}
        }
        //Find total_average
        for( int i =0; i < 5; i++){
          MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // print results
        if ( my_rank == 0) {
          output(n_spins, mcs, temperature, total_average);
        }
    }

    //    cout << "analytical Eavg : " << -32.*sinh(8.)/(12. + 4.*cosh(8.))/4. << endl;
    //    cout << "analytical Evar : " << 1024.*(1.+ 3.*cosh(8.)) / pow(12. + 4*cosh(8.),2)/4. << endl;
    //    cout << "analytical Mavg : " << 0. << endl;
    //    cout << "analytical Mvar : " << 32.*(exp(8.)+1.)/(12.+4.*cosh(8.))/4. - pow(8.*(exp(8.)+2.)/(12.+ 4*cosh(8.)) ,2)/4 << endl;
    //    cout << "analytical Mabs average : " << 8.*(exp(8.)+2.)/(12.+ 4*cosh(8.))/4. << endl;
    //cout << "Accepted moves: "<< acceptedmoves <<endl;
    //cout << "Total moves: "<<mcs*n_spins*n_spins << endl;
    //cout << "Percent accepted: "<<(double)acceptedmoves/(double)(mcs*n_spins*n_spins)*100. << endl;
    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file
    //End MPI
    MPI_Finalize ();
    return 0;
}


// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
                double& final_temp, double& temp_step)
{
    cout << "Number of Monte Carlo trials =";
    cin >> mcs;
    cout << "Lattice size or number of spins (x and y equal) =";
    cin >> n_spins;
    cout << "Initial temperature with dimension energy=";
    cin >> initial_temp;
    cout << "Final temperature with dimension energy=";
    cin >> final_temp;
    cout << "Temperature step with dimension energy=";
    cin >> temp_step;
} // end of function read_input


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix,
                double& E, double& M)
{
    long idum;
    idum = -1;
    double randomnumber;
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            if (temperature <1.5) spin_matrix[y][x] = 1; // spin orientation for the ground state
            //spin_matrix[y][x] = 1; // spin orientation for the ground state
            randomnumber = ran0(&idum);
            if(randomnumber >= 0.5) spin_matrix[y][x] = 1;
            if(randomnumber < 0.5) spin_matrix[y][x] = -1;
            M +=  (double) spin_matrix[y][x];
        }
    }
    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
                    (spin_matrix[periodic(y,n_spins,-1)][x] +
                    spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
    // loop over all spins
    //accepted = 0;
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            int ix = (int) (ran1(&idum)*(double)n_spins);
            int iy = (int) (ran1(&idum)*(double)n_spins);
            int deltaE =  2*spin_matrix[iy][ix]*
                    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                    spin_matrix[periodic(iy,n_spins,-1)][ix] +
                    spin_matrix[iy][periodic(ix,n_spins,1)] +
                    spin_matrix[periodic(iy,n_spins,1)][ix]);
            if ( ran1(&idum) <= w[deltaE+8] ) {
                spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
                //accepted +=1;
            }
        }
    }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *average)
{
    double norm = 1/((double) (mcs));  // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;
} // end output function
