#include "mpi_grid.h"

MPI_Grid::MPI_Grid(int myId0, int nProc0, int nx0, int ny0)
{
    myId = myId0;
    nProc = nProc0;
    nx = nx0;
    ny = ny0;

    mpi_decompose();
}

/**
 * Decomposes the domain for MPI
 */
void MPI_Grid::mpi_decompose(void){
    double scale = ((float)nx)/ny;
    double sDiff, s0;
    nRow=1, nCol=nProc;

    sDiff = std::abs(1.0/nProc - scale);
    for(int i=2; i<=nProc; i++){
        if(nProc%(i) == 0){
            s0 = std::abs(i/float(nProc/i) - scale);
            if(s0 <= sDiff){
                sDiff = s0;
                nRow = i;
                nCol = nProc/i;
            }
        }
    }

    // MPI grid indexes
    rowIdx = int(myId/nCol);
    colIdx = int(myId%nCol);

    //Neighbor MPI Ids
    xpId = rowIdx*nCol + (colIdx+1)%(nCol);
    xmId = rowIdx*nCol + (nCol+colIdx-1)%(nCol);
    // The y indexes may appear swapped but this is
    // because the origin of the grid is upper left
    ypId = nCol*((rowIdx+1)%nRow) + colIdx;
    ymId = nCol*((nRow+rowIdx-1)%nRow) + colIdx;
    
    // Local domain size, recall x is determined by columns and y by rows
    lx = int((colIdx+1.0)/nCol*nx) - int(float(colIdx)/nCol*nx);
    ly = int((rowIdx+1.0)/nRow*ny) - int(float(rowIdx)/nRow*ny);
}

/**
 * Gets the local domain size of the MPI process
 * @param lx0, ly0 integer pointers for local domain node size
 */
void MPI_Grid::get_domain_size(int *lx0, int *ly0){
    *lx0 = lx;
    *ly0 = ly;
}
/**
 * Accessor global domain size
 * @param nx0, ny0 integer pointers for local domain node size
 */
void MPI_Grid::get_global_domain(int *nx0, int *ny0){
    *nx0 = nx;
    *ny0 = ny;
}

/**
 * Method to determine which sides of local domain
 * are interior within the global domain
 * @return boolean array of length 4, with true if interior side
 * [top, right, bottom, left]. Origin is in bottom left.
 */
void MPI_Grid::get_interior_bounds(bool* interiorSides){
    for(int i=0; i<4; i++)
        interiorSides[i] = true;

    // Check if local domain is on global domain boundary
    if(rowIdx == nRow-1){
        interiorSides[0] = false;
    }
    if(colIdx == nCol-1){
        interiorSides[1] = false;
    }
    if(rowIdx == 0){
        interiorSides[2] = false;
    }
    if(colIdx == 0){
        interiorSides[3] = false;
    }
}
/**
 * Gets the location of the local MPI domain with respect to the
 * the global domain.
 * @param mpiIdx - integer index of process we would like the get the domain of
 * @param lOrigin - integer array of domains local origin [nc, nr]
 * @param lSize - integer array of mpi tasks local matrix size [ly, lx]
 */
void MPI_Grid::get_proc_domain(int mpiIdx, int *lOrigin, int *lSize){

    // MPI grid indexes
    int rowIdx0 = int(mpiIdx/nCol);
    int colIdx0 = int(mpiIdx%nCol);

    // Local domain size, recall x is determined by columns and y by rows
    lSize[1] = int((colIdx0+1.0)/nCol*nx) - int(float(colIdx0)/nCol*nx);
    lSize[0] = int((rowIdx0+1.0)/nRow*ny) - int(float(rowIdx0)/nRow*ny);

    lOrigin[1] = 0;
    for(int i=0; i<colIdx0; i++){
        lOrigin[1] += int((i+1.0)/nCol*nx) - int(float(i)/nCol*nx);
    }
    lOrigin[0] = 0;
    for(int i=0; i<rowIdx0; i++){
        lOrigin[0] += int((i+1.0)/nRow*ny) - int(float(i)/nRow*ny);
    }
}

/**
 * Accessores for MPI neighbor indexes
 */
int MPI_Grid::get_yp_id(void){
    return ypId;
}
int MPI_Grid::get_ym_id(void){
    return ymId;
}
int MPI_Grid::get_xp_id(void){
    return xpId;
}
int MPI_Grid::get_xm_id(void){
    return xmId;
}
/**
 * Accessor for MPI id
 */
int MPI_Grid::get_id(void){
    return myId;
}

int MPI_Grid::get_nproc(void){
    return nProc;
}