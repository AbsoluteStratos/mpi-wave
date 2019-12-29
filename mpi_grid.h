#include "common.h"

#if !defined(_MPI_GRID_H)
#define _MPI_GRID_H

#if defined(__MPI__)
#include "mpi.h"
#endif /*if defined(__MPI__)*/

class MPI_Grid{
	int nProc;
    int myId;

    // Number of nodes
    int nx, lx;
    int ny, ly;

    // Local MPI info
    int nRow, nCol;
    int rowIdx, colIdx;
    int xpId, xmId, ypId, ymId;

public:
	MPI_Grid(){ // default constructor
	    nProc = 1;
	    myId = 0;
        
        nx = 100;
        ny = 100;

        mpi_decompose();
	}
    MPI_Grid(int, int, int, int);

    void get_domain_size(int*, int*);
    void get_global_domain(int*, int*);
    void get_interior_bounds(bool*);
    void get_proc_domain(int, int*, int*);
    int get_yp_id(void);
    int get_ym_id(void);
    int get_xp_id(void);
    int get_xm_id(void);
    int get_id(void);
    int get_nproc(void);
private:
    void mpi_decompose(void);
};

#endif 