/*
* C++ solver for the 2D wave equation through hetergenous media  
* with MPI parallelization.
* Primary author: Nicholas Geneva (ngeneva at nd.edu)
*
* To compile in parallel: 
* mpicxx *cpp -std=c++11 -lboost_system -lboost_filesystem  -D__MPI__
*/
#include "common.h"
#include "mpi_grid.h"
#include "boundary.h"
#include "data_io.h"
#include "source.h"
#include "finite_difference.h"

#if defined(__MPI__)
#include  "mpi.h"
#endif /*if defined(__MPI__)*/

using namespace std;

int main(int argc, char *argv[])
{
      int myId, nProc;
      // ===== MODEL PARAMETERS ======= //
      int nx = 250, ny = 250;
      // 0 = periodic, 1 = dirichlet, 2 = open
      int bType = 0;

      /* Time-step variables */
      // double tStart= 0.0;
      double tEnd = 20;
      double dt = 0.01;

      int saveItr = 25;

      /* Timing variables */
      double starttime = 0.0;
      double endtime = 0.0;

#if defined(__MPI__)
      MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &nProc);
      MPI_Comm_rank(MPI_COMM_WORLD, &myId);

      MPI_Grid mpi_grid(myId, nProc, nx, ny);
      starttime = MPI_Wtime();
#else
      myId = 0;
      nProc = 1;
      MPI_Grid mpi_grid(myId, nProc, nx, ny);
#endif

      DataIO dataio(mpi_grid);
      Source source(mpi_grid);
      int sX = nx/4; int sY = ny/4; 
      double sourceMag = 0.1/(dt*dt)*(nx/50.)*(ny/50.);
      // Point source: xpos, ypos, mag, duration
      source.set_point_source(sX, sY, sourceMag, 2.0);
      // Harmonic source: xpos, ypos, mag, duration, freq.
      // source.set_harmonic_source(sX, sY, sourceMag, 5.0, 2.0);

      Finite_Difference fd(mpi_grid, 1.0/nx, 1.0/ny, dt);

      int bnds[4];
      if(bType == 0){
            // Periodic
            bnds[0]=0; bnds[1]=0; bnds[2]=0; bnds[3]=0;
      }else if(bType == 1){
            // Dirichlet
            bnds[0]=1; bnds[1]=1; bnds[2]=1; bnds[3]=1;
      }else if(bType == 2){
            // Open
            bnds[0]=2; bnds[1]=2; bnds[2]=2; bnds[3]=2;
      }else{
            cout << "Boundary type not supported. Defaulting to periodic..." << endl;
            bnds[0]=0; bnds[1]=0; bnds[2]=0; bnds[3]=0;
      }
      Boundary boundary(mpi_grid, bnds);
      boundary.set_discretization(1.0/nx, 1.0/ny, dt);

      int lx, ly;
      mpi_grid.get_domain_size(&lx, &ly);
      cout << myId << " Domain size: " << lx <<"," << ly<<endl;

      // Start wall-clock execution timer
      starttime = MPI_Wtime();

      // Allocate state matrix
      double** U;
      alloc_2d_array(ly+2, lx+2, U);
      for(int i=1; i<=ly; i++){
            for(int j=1; j<=lx; j++){
                  U[i][j] = 0;
            }  
      }
      // Allocate previous time-step state matrix
      double** U0;
      alloc_2d_array(ly, lx, U0);
      for(int i=0; i<ly; i++){
            for(int j=0; j<lx; j++){
                  U0[i][j] = 0;
            }  
      }
      // Allocate source field
      double** S;
      alloc_2d_array(ly, lx, S);

      // Allocate and read material field
      double** K;
      alloc_2d_array(ly, lx, K);
      dataio.load_material(K, "perm_gen/mat.bin");
      // Uniform
      // for(int i=0; i<ly; i++){
      //       for(int j=0; j<lx; j++){
      //             K[i][j] = 0.1;
      //       }  
      // }
      boundary.set_material_field(K);

      // Block to make sure all processes get thier stuff allocated
      MPI_Barrier(MPI_COMM_WORLD);
      
      // ====================== Main Time Loop ====================== //
      double currT = 0;
      for(int t=0; t<int(tEnd/dt); t++){

            // First deal with boundaries
            boundary.set_boundaries(U);

            // Update source field
            source.get_source(S, currT);

            // Data log
            if(t%saveItr == 0){
                  if(myId == 0){
                        cout << "Saving current time-step: "<< round(currT*100)/100 << endl;
                  }     
                  dataio.save_state(U, currT);
            }

            // Forward euler
            fd.forward_euler(U, U0, K, S);
            currT+=dt;
      }

      MPI_Barrier(MPI_COMM_WORLD);
      endtime = MPI_Wtime();

      if(myId == 0){
            cout << myId << ": Took " << round(1000*(endtime-starttime))/1000 << " sec. to execute." << endl;
      }
            
#if defined(__MPI__)
      MPI_Finalize();
#endif     
}

