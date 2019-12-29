#include "finite_difference.h"

/**
 * Constructor
 */
 Finite_Difference::Finite_Difference(MPI_Grid mpi_grid0, double dx0, double dy0, double dt0){
    mpi_grid = mpi_grid0;
    // Discretization parameters
    dx = dx0;
    dy = dy0;
    dt = dt0;

    // Allocate array for temporary storage of the next time-step
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    alloc_2d_array(ly, lx, Udt);
 }

/**
 * Forward time-integration method, predicted time-step t+dt.
 * Ghost nodes should already be updated before calling this!
 * @param U - 2D array of state variables [ly+2 x lx+2]
 * @param U0 - 2D array of state variables at t-dt [ly x lx]
 * @param K - 2D array of material property [ly x lx]
 * @param S - 2D array of source field [ly x lx]
 */
void Finite_Difference::forward_euler(double**& U, double**& U0, double**& K, double**& S){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    double c0y = dt/dy;
    double c0x = dt/dx;
    
    // Finite diference here
    // See report for details
    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            // Note: when accessing U, we must account for 1 layer of ghost nodes
            Udt[i][j] = pow(K[i][j]*c0y, 2)*U[i+2][j+1] + pow(K[i][j]*c0x, 2)*U[i+1][j+2] + 
                2.0*(1.0-pow(K[i][j]*c0x, 2) - pow(K[i][j]*c0y, 2))*U[i+1][j+1] + pow(K[i][j]*dt, 2)*S[i][j] + 
                pow(K[i][j]*dt/dy, 2)*U[i][j+1] + pow(K[i][j]*c0x, 2)*U[i+1][j]-U0[i][j];
        }
    }
    // Update past time-step array
    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            U0[i][j] = U[i+1][j+1];
        }
    }
    // Update state array
    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            U[i+1][j+1] = Udt[i][j];
        }
    }
        
}

/**
 * Modifiers for discreatization parameters
 */
void Finite_Difference::set_dx(double dx0){
    dx = dx0;
}
void Finite_Difference::set_dy(double dy0){
    dy = dy0;
}
void Finite_Difference::set_dt(double dt0){
    dt = dt0;
}