#if !defined(_SOURCE_H)
#define _SOURCE_H

#include "common.h"
#include "mpi_grid.h"

#if defined(__MPI__)
#include "mpi.h"
#endif /*if defined(__MPI__)*/

class Source{

    MPI_Grid mpi_grid;
    MPI_Status  dstatus;
    int sourceNx, sourceNy;
    double sourceMaxT, sourceMag, sourcePeriod;

public:
    Source(){ // default constructor
        //TODO
	}

    Source(MPI_Grid);

    void get_source(double**&, double);
    void set_no_source(void);
    void set_point_source(int, int, double, double);
    void set_harmonic_source(int, int, double, double, double);

private:

    // Source function pointer
    void (Source::*sourceFunc)(double**&, double);

    void zero_source(double**&, double);
    void point_source(double**&, double);
    void harmonic_source(double**&, double);
};

#endif 