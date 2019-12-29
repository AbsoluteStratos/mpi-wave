#if !defined(_DATA_IO_H)
#define _DATA_IO_H

#include "common.h"
#include "mpi_grid.h"
#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>

#if defined(__MPI__)
#include "mpi.h"
#endif /*if defined(__MPI__)*/

class DataIO{

    MPI_Grid mpi_grid;
    MPI_Status  dstatus;
    std::ofstream ofp;
    std::ifstream ifp;

    int nProc;
    int nx, ny;

public:
    DataIO(){ // default constructor
        //TODO
        nProc = 1;
	}

    DataIO(MPI_Grid);

    void save_state(double**&, double);
    void load_material(double**&, const char*);
    // void save_mpi_decomp(double);

private:
    void save_array_binary(double**&, int, int, const char*);
    void load_array_binary(double**&, int, int, const char*);
    void save_array_text(double**&,  int, int, const char*);

};

#endif 