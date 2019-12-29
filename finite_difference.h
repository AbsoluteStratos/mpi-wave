#if !defined(_FINITE_DIFF_H)
#define _FINITE_DIFF_H

#include "common.h"
#include "mpi_grid.h"

#if defined(__MPI__)
#include "mpi.h"
#endif /*if defined(__MPI__)*/

class Finite_Difference{

    MPI_Grid mpi_grid;
    int nx, ny;
    double dx, dy, dt;
    double** Udt;
public:
    Finite_Difference(){ // default constructor
        //TODO
	}

    Finite_Difference(MPI_Grid, double, double, double);

    void forward_euler(double**&, double**&, double**&, double**&);

    // Modifiers incase adaptive time-step size/ mesh
    void set_dx(double);
    void set_dy(double);
    void set_dt(double);
};

#endif