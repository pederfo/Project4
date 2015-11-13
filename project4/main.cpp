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
void initialize(int, double, int **, double&, double&, bool ordered);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *, long &);
// prints to file the results of the calculations
void output(int, int, double, double *);

int main(int argc, char* argv[])
{
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs, my_rank, numprocs,method;
    double w[17],initial_temp, final_temp, E, M, temp_step;
    bool ordered;
    int counter= 0;


    // MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Read in output file, abort if there are too few command-line arguments
    if (my_rank == 0 && argc <= 1) {
        cout << "Bad Usage: " << argv[0] <<
                " read output file" << endl;
        exit(1);
    }
    if (my_rank == 0 && argc > 1) {
        outfilename=argv[1];
        ofile.open(outfilename);
    }
    /*
    Define parameters for this run of the program. Note the method parameter. This paramter
    defines what information the user wants the program to print, either to file, to terminal, or both.
    Descriptions of the methods can be found in the readme file located in the same folder as the program.
    Note also the "ordered" parameter defining whether the user wants an ordered system at the beginning
    or an unordered one. This is controlled simply with either true or false statements
    */

    n_spins = 80; mcs = 1000000;  initial_temp = 2.0; final_temp = 2.7; temp_step =0.01; method = 2,ordered=true;

    // Broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //Allocate memory for spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));


    //Array for energies.This can only be used for one temperature at a time
    double energies[mcs];
    for(int i=0;i<mcs;i++) energies[i] = 0;

    // Determine the values for all temperature steps and assign each to their respective processors.
    bool extra = false;
    int num_temps = (final_temp - initial_temp)/temp_step+ 1;
    int my_n = num_temps/numprocs;
    int rest = num_temps%numprocs;
    int startindex = 0;
    if(my_rank < rest){
        my_n+=1;
        extra = true;
    }
    if(extra){
        startindex = my_rank*my_n;
    }
    else{
        startindex = num_temps -my_n*(numprocs - my_rank);
    }
    //Define the average variables we will be filling
    double average[num_temps][5];
    //Assign arrays for counting the number of accepted moves for each temperature step.
    long accepted[num_temps];
    for( int i = 0; i < num_temps; i++) accepted[i] = 0.;
    long acceptedmoves[num_temps];
    for( int i = 0; i < num_temps; i++) acceptedmoves[i] = 0.;
    long acceptedmoves_global[num_temps];
    for( int i = 0; i <= num_temps; i++) acceptedmoves_global[i] = 0.;

    //Fill the temperature array with all the temperatures for this run.
    double temperature[num_temps];
    for(int i=0;i<num_temps;i++){
        temperature[i] = initial_temp +i*temp_step;
    }

    //Brute force method for when to start saving expectation values.
    int startcount = 300000;

    //double energies[mcs];
    idum = -1-my_rank; // random starting point
    for ( int i = startindex; i< startindex+my_n; i++){
        cout << "Running calculation for temperature = " << temperature[i] << endl;
        //Initialise energy and magnetization
        E = M = 0.;
        //Setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature[i]);
        //Initialise array for expectation values
        for( int j = 0; j < 5; j++) average[i][j] = 0.;
        initialize(n_spins, temperature[i], spin_matrix, E, M,ordered);

        // start Monte Carlo computation
        for (int cycles = 0; cycles <= mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w,accepted[i]);
            // update expectation values
            average[i][0] += E;    average[i][1] += E*E;
            average[i][2] += M;    average[i][3] += M*M; average[i][4] += fabs(M);
            if(cycles >= startcount){
                energies[cycles] = E;
                if((my_rank==0) & (method==4)){
                    ofile << setw(15) << energies[cycles] << endl;
                }
                acceptedmoves[i] += accepted[i];
                accepted[i] = 0;

            }

            if((my_rank==0)&(method==3)& (counter ==1000)){
                ofile << setw(15) << setprecision(8) << cycles;
                output(n_spins, cycles, temperature[i], average[i]);
                counter = 0;
            }

            counter +=1;
            accepted[i] = 0;
        }

    }

    if((my_rank == 0) & (method==0)){
        cout << "analytical Eavg : " << -32.*sinh(8.)/(12. + 4.*cosh(8.))/4. << endl;
        cout << "analytical Evar : " << 1024.*(1.+ 3.*cosh(8.)) / pow(12. + 4*cosh(8.),2)/4. << endl;
        cout << "analytical Mavg : " << 0. << endl;
        cout << "analytical Mvar : " << 32.*(exp(8.)+1.)/(12.+4.*cosh(8.))/4. - pow(8.*(exp(8.)+2.)/(12.+ 4*cosh(8.)) ,2)/4 << endl;
        cout << "analytical Mabs average : " << 8.*(exp(8.)+2.)/(12.+ 4*cosh(8.))/4. << endl;
        output(n_spins, mcs, temperature[0], average[0]);
    }

    //Sort all acceptemoves from all processors into one global array
    MPI_Allreduce(acceptedmoves,acceptedmoves_global,num_temps,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

    //Sort all average values into a global array
    double global_average[num_temps][5];
    for(int i = 0; i<num_temps;i++){
        for(int j=0;j<5;j++){
            global_average[i][j] = 0;
        }
    }
    double temperaturelist[num_temps];
    for(int i=0;i<num_temps;i++){
           temperaturelist[i] = temperature[i];
    }

    for(int i=0;i<num_temps;i++){
        MPI_Allreduce(average[i],global_average[i],num_temps,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    for(int i=0;i<num_temps;i++){
        cout << temperaturelist[i]<< endl;
    }
    /*This method prints the temperature step and the numebr of accepted moves for that step.
    The user must note the number of cycles and the number of spins themselves if they want
    to plot the percentage of the total moves. (For reference the total numer of moves=n_spins*n_spins*mcs)
    */
    if((my_rank==0) & (method==1)){
        for(int i=0;i<num_temps;i++){
            ofile << setw(15) << setprecision(8) << initial_temp + i*temp_step <<"\t"<< acceptedmoves_global[i] << endl;
        }
    }

    /*This method prints the expectation values
    for every temperature in the interval set by the user*/
    if((method==2) & (my_rank==0)){
        for(int i=0;i<num_temps;i++){
            output(n_spins,mcs,temperaturelist[i],global_average[i]);
        }
    }


    //cout << "Accepted moves: "<< acceptedmoves <<endl;
    //cout << "Total moves: "<<mcs*n_spins*n_spins << endl;
    //cout << "Percent accepted: "<<(double)acceptedmoves/(double)(mcs*n_spins*n_spins)*100. << endl;

    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file
    //End MPI
    MPI_Finalize ();
    return 0;
}
//End of main program



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
                double& E, double& M,bool ordered)
{
    long idum;
    idum = -1;
    double randomnumber;
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            //if (temperature <1.5) spin_matrix[y][x] = 1; // Sets state ordered if temperature is very low

            // Ordered state
            if(ordered==true) spin_matrix[y][x] = 1;
            //Disordered state
            if(ordered==false){
                randomnumber = ran0(&idum);
                if(randomnumber >= 0.5) spin_matrix[y][x] = 1;
                if(randomnumber < 0.5) spin_matrix[y][x] = -1;
            }
            //Fill array
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


//Function to do the metropolis sampling
void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, long &accepted)
{
    // loop over all spins
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
                accepted +=1;
            }
        }
    }
} // end of Metropolis sampling over spins

//Function that prints to file
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
