#include "source.h"

Source::Source(MPI_Grid mpi_grid0){
    mpi_grid = mpi_grid0;
    // Default source function is zero source
    sourceFunc = &Source::zero_source;
}

/**
 * Gets source field for current time t. Public wrapper
 * for different source functions
 * @param S - 2D source field array for the local MPI domain
 * @param t - current time
 */
void Source::get_source(double**& S, double t){
    (this->*sourceFunc)(S, t);
}

/**
 * Sets source field to zero over entire domain (default)
 */
void Source::set_no_source(void){
    sourceFunc = &Source::zero_source;
}

/**
 * Sets source function to point source
 * @param xLoc, yLoc  - integers of sources node location on 
 *      the GLOBAL domain with size [0, nx-1]x[0,ny-1].
 * @param mag - magnitude of point source
 * @param tMax - max time the source should be active
 */
void Source::set_point_source(int xLoc, int yLoc, double mag, double tMax){
    int lSize[2];
    int lOrigin[2];
    int myId = mpi_grid.get_id();
    // Remember orgin and size are given [col, row]
    mpi_grid.get_proc_domain(myId, lOrigin, lSize);

    // We first need to check if the source is in our local MPI tasks domain
    if(yLoc >= lOrigin[0] and yLoc-lOrigin[0] < lSize[0]){
        if(xLoc >= lOrigin[1] and xLoc-lOrigin[1] < lSize[1]){
            // Set local location of source
            sourceNx = xLoc-lOrigin[1];
            sourceNy = yLoc-lOrigin[0];

            sourceMag = mag;
            sourceMaxT = tMax;

            sourceFunc = &Source::point_source;
        }else{
            // If not in domain just set it to zero
            sourceFunc = &Source::zero_source;
        }
    }else{
        // If not in domain just set it to zero
        sourceFunc = &Source::zero_source;
    }
}

/**
 * Sets source function to point source
 * @param xLoc, yLoc  - integers of sources node location on
 *      the GLOBAL domain with size [0, nx-1]x[0,ny-1].
 * @param mag - magnitude of point source
 * @param tMax - max time the source should be active
 * @param tPeriod - period of harmonic source w.r.t. time
 */
void Source::set_harmonic_source(int xLoc, int yLoc, double mag, double tMax, double tPeriod){
    int lSize[2];
    int lOrigin[2];
    int myId = mpi_grid.get_id();
    // Remember orgin and size are given [col, row]
    mpi_grid.get_proc_domain(myId, lOrigin, lSize);

    // We first need to check if the source is in our local MPI tasks domain
    if(yLoc >= lOrigin[0] and yLoc-lOrigin[0] < lSize[0]){
        if(xLoc >= lOrigin[1] and xLoc-lOrigin[1] < lSize[1]){
            // Set local location of source
            sourceNx = xLoc-lOrigin[1];
            sourceNy = yLoc-lOrigin[0];

            sourceMag = mag;
            sourceMaxT = tMax;
            sourcePeriod = tPeriod;

            sourceFunc = &Source::harmonic_source;
        }else{
            // If not in domain just set it to zero
            sourceFunc = &Source::zero_source;
        }
    }else{
        // If not in domain just set it to zero
        sourceFunc = &Source::zero_source;
    }
}

/**
 * Zero source function
 * @param S - 2D source field array for the local MPI domain
 * @param t - current time
 */
void Source::zero_source(double**& S, double t){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            S[i][j] = 0;
        }  
    }
}

/**
 * Point source function
 * @param S - 2D source field array for the local MPI domain
 * @param t - current time
 */
void Source::point_source(double**& S, double t){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            S[i][j] = 0;
        }  
    }

    // If time is still less than max source time
    if(t <= sourceMaxT){
        S[sourceNy][sourceNx] = sourceMag;
    }
}

/**
 * Harmonic point source function
 * @param S - 2D source field array for the local MPI domain
 * @param t - current time
 */
void Source::harmonic_source(double**& S, double t){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    for(int i=0; i<ly; i++){
        for(int j=0; j<lx; j++){
            S[i][j] = 0;
        }  
    }

    // If time is still less than max source time
    if(t <= sourceMaxT){
        S[sourceNy][sourceNx] = sourceMag*cos(2*PI*(t/sourcePeriod));
    }
}