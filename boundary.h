#include "common.h"
#include "mpi_grid.h"

#if defined(__MPI__)
#include "mpi.h"
#endif /*if defined(__MPI__)*/

class Boundary{

    MPI_Grid mpi_grid;
    MPI_Status  bstatus;
    MPI_Request brequestYp, brequestXp, brequestYm, brequestXm;
    bool mpiComms[4];
    double bndValues[4];
    double dt, dx, dy;
    double **K;
    double *ypSendBuff, *xpSendBuff, *ymSendBuff, *xmSendBuff;
    double *ypRecvBuff, *xpRecvBuff, *ymRecvBuff, *xmRecvBuff;

public:
    Boundary(){ // default constructor
        //TODO
	}

    Boundary(MPI_Grid, int*);

    void set_discretization(double,double,double);
    void set_material_field(double**&);
    void set_boundaries(double**&);
    void mpi_wait_boundaries(void);
    
private:
    void init_boundaries(int*);
    void mpi_send_boundaries(double**&);
    
    // Boundary function pointers
    void (Boundary::*bxp)(double**&);
    void (Boundary::*bxm)(double**&);
    void (Boundary::*byp)(double**&);
    void (Boundary::*bym)(double**&);

    // Boundary functions
    void interiorXp(double**&);
    void interiorXm(double**&);
    void interiorYp(double**&);
    void interiorYm(double**&);

    void dirichletXp(double**&);
    void dirichletXm(double**&);
    void dirichletYp(double**&);
    void dirichletYm(double**&);

    void localPeriodicXp(double**&);
    void localPeriodicXm(double**&);
    void localPeriodicYp(double**&);
    void localPeriodicYm(double**&);

    void openXp(double**&);
    void openXm(double**&);
    void openYp(double**&);
    void openYm(double**&);
};