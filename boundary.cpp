#include "boundary.h"

/**
 * Constructor
 */
Boundary::Boundary(MPI_Grid mpi_grid0, int* bndryArr){
    mpi_grid = mpi_grid0;

    for(int i=0; i<4; i++){
        mpiComms[i] = false;
        bndValues[i] = 0;
    }

    brequestYp = MPI_REQUEST_NULL;
    brequestXp = MPI_REQUEST_NULL;
    brequestYm = MPI_REQUEST_NULL;
    brequestXm = MPI_REQUEST_NULL;

    init_boundaries(bndryArr);
}

/**
 * Initializes the boundary types by setting the boundary
 * pointers to the correct function. If MPI is needed
 * this will also allocate necessary buffers.
 * @param bndryArr - int array of 4 values of the boundary
 *          type of the global domain. Index representation:
 *          0) Periodic, 1) Dirichlet, 2) Open.
 *          Order: [y+, x+, y-, x-]
 */
void Boundary::init_boundaries(int* bndryArr){
    bool intrBnds[4];
    int lx, ly;
    int myId = mpi_grid.get_id();
    mpi_grid.get_domain_size(&lx, &ly);

    // First get a list of interious boundarys of local domain
    mpi_grid.get_interior_bounds(intrBnds);
    
    // Top boundary
    if(bndryArr[0] == 0 || intrBnds[0] == true){ // Periodic/Interior
        if(myId == mpi_grid.get_yp_id()){
            byp = &Boundary::localPeriodicYp;
        }else{
            byp = &Boundary::interiorYp;
            mpiComms[0] = true;
            alloc_1d_array(lx, ypSendBuff);
            alloc_1d_array(lx, ypRecvBuff);
        }
    }else if(bndryArr[0] == 1){ // Dirichlet
        byp = &Boundary::dirichletYp;
    }else if(bndryArr[0] == 2){ // Open
        byp = &Boundary::openYp;
    }

    // Right boundary
    if(bndryArr[1] == 0 || intrBnds[1] == true){ // Periodic/Interior
        if(myId == mpi_grid.get_xp_id()){
            bxp = &Boundary::localPeriodicXp;
        }else{
            bxp = &Boundary::interiorXp;
            mpiComms[1] = true;
            alloc_1d_array(ly, xpSendBuff);
            alloc_1d_array(ly, xpRecvBuff);
        }
    }else if(bndryArr[1] == 1){ // Dirichlet
        bxp = &Boundary::dirichletXp;
    }else if(bndryArr[1] == 2){ // Open
        bxp = &Boundary::openXp;
    }

    // Lower boundary
    if(bndryArr[2] == 0 || intrBnds[2] == true){ // Periodic/Interior
        if(myId == mpi_grid.get_ym_id()){
            bym = &Boundary::localPeriodicYm;
        }else{
            bym = &Boundary::interiorYm;
            mpiComms[2] = true;
            alloc_1d_array(lx, ymSendBuff);
            alloc_1d_array(lx, ymRecvBuff);
        }
    }else if(bndryArr[2] == 1){ // Dirichlet
        bym = &Boundary::dirichletYm;
    }else if(bndryArr[2] == 2){ // Open
        bym = &Boundary::openYm;
    }

    // Left boundary
    if(bndryArr[3] == 0 || intrBnds[3] == true){ //Interior Boundary
        if(myId == mpi_grid.get_xm_id()){
            bxm = &Boundary::localPeriodicXm;
        }else{
            bxm = &Boundary::interiorXm;
            mpiComms[3] = true;
            alloc_1d_array(ly, xmSendBuff);
            alloc_1d_array(ly, xmRecvBuff);
        }
    }else if(bndryArr[3] == 1){ // Dirichlet
        bxm = &Boundary::dirichletXm;
    }else if(bndryArr[3] == 2){ // Open
        bxm = &Boundary::openXm;
    }
}

/**
 * Used to set discretization parameters of the domain. Should be called
 * when using the open boundary condition.
 * @param dx0 - double of mesh size in the x direction
 * @param dy0 - double of mesh size in the y direction
 * @param dt0 - double of time-step size
 */
void Boundary::set_discretization(double dx0, double dy0, double dt0){
    dx = dx0;
    dy = dy0;
    dt = dt0;
}

/**
 * Used to set material array pointer. Should be called
 * when using the open boundary condition.
 * @param K0 - 2D array of material property for local MPI domain
 */
void Boundary::set_material_field(double**& K0){
    K = K0;
}

void Boundary::set_boundaries(double**& A){
    // First send all data to neighors, async
    #if defined(__MPI__)
        mpi_send_boundaries(A);
    #endif

    // Call boundary function pointers
    (this->*byp)(A);
    (this->*bxp)(A);
    (this->*bym)(A);
    (this->*bxm)(A);
}

/**
 * Sends needed MPI boundaries to the respective neighbors
 * A local node sends data if the boundary is interior or
 * there is a periodic global boundary.
 * @param A - 2D tensor of state variables in local domain
 */
void Boundary::mpi_send_boundaries(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    // Before sending make sure all past messages
    // have been properly recieved.
    mpi_wait_boundaries();

    // Send yp boundary
    if(mpiComms[0]==true){
        for(int i=0; i<lx; i++){
            ypSendBuff[i] = A[ly][i+1];
        }
        int dest = mpi_grid.get_yp_id();
        // std::cout << mpi_grid.get_id() << ": Yp Sending data to: "<< dest << std::endl;
        MPI_Isend(&ypSendBuff[0], lx, MPI_DOUBLE, dest, 1231, MPI_COMM_WORLD, &brequestYp);
    }
    // Send xp boundary
    if(mpiComms[1]==true){
        for(int j=0; j<ly; j++){
            xpSendBuff[j] = A[j+1][lx];
        }
        int dest2 = mpi_grid.get_xp_id();
        // std::cout << mpi_grid.get_id() << ": Xp Sending data to: "<< dest2 << std::endl;
        MPI_Isend(&xpSendBuff[0], ly, MPI_DOUBLE, dest2, 1232, MPI_COMM_WORLD, &brequestXp);
    }
    // Send ym boundary
    if(mpiComms[2]==true){
        for(int i=0; i<lx; i++){
            ymSendBuff[i] = A[1][i+1];
        }
        int dest3 = mpi_grid.get_ym_id();
        // std::cout << mpi_grid.get_id() << ": Ym Sending data to: "<< dest3 << std::endl;
        MPI_Isend(&ymSendBuff[0], lx, MPI_DOUBLE, dest3, 1233, MPI_COMM_WORLD, &brequestYm);
    }
    // Send xm boundary
    if(mpiComms[3]==true){
        for(int j=0; j<ly; j++){
            xmSendBuff[j] = A[j+1][1];
        }
        int dest4 = mpi_grid.get_xm_id();
        // std::cout << mpi_grid.get_id() << ": Xm Sending data to: "<< dest4 << std::endl;
        MPI_Isend(&xmSendBuff[0], ly, MPI_DOUBLE, dest4, 1234, MPI_COMM_WORLD, &brequestXm);
    }
}

/**
 * Sends needed MPI boundaries to the respective neighbors
 * A local node sends data if the boundary is interior or
 * there is a periodic global boundary.
 * @param 2D tensor of state variables in local domain
 */
void Boundary::mpi_wait_boundaries(void){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

     // Send yp boundary
    if(mpiComms[0]==true){
        MPI_Wait(&brequestYp, &bstatus);
    }
    // Send xp boundary
    if(mpiComms[1]==true){
        MPI_Wait(&brequestXp, &bstatus);
    }
    // Send ym boundary
    if(mpiComms[2]==true){
        MPI_Wait(&brequestYm, &bstatus);
    }
    // Send xm boundary
    if(mpiComms[3]==true){
        MPI_Wait(&brequestXm, &bstatus);
    }
}

/**
 * MPI Interior (also used for periodic) boundaries
 * Used to recieve ghost nodes from neighboring process
 * and set them in the local state variable array.
 */
void Boundary::interiorYp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    int dest = mpi_grid.get_yp_id();
    // Remember tag should be for a Ym send
    // std::cout << mpi_grid.get_id() << ": Yp Waiting data from: "<< dest << std::endl;
    MPI_Recv(ypRecvBuff, lx, MPI_DOUBLE, dest, 1233, MPI_COMM_WORLD, &bstatus);
    // std::cout << mpi_grid.get_id() << ": Recieved data from: "<< dest << std::endl;
    // Set A's ghost nodes
    for(int i=0; i<lx; i++){
        A[ly+1][i+1] = ypRecvBuff[i];
    }
}
void Boundary::interiorXp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    int dest = mpi_grid.get_xp_id();
    // Remember tag should be for a Xm send
    // std::cout << mpi_grid.get_id() << ": Xp Waiting data from: "<< dest << std::endl;
    MPI_Recv(xpRecvBuff, ly, MPI_DOUBLE, dest, 1234, MPI_COMM_WORLD, &bstatus);
    // std::cout << mpi_grid.get_id() << ": Recieved data from: "<< dest << std::endl;
    // Set A's ghost nodes
    for(int j=0; j<ly; j++){
        A[j+1][lx+1] = xpRecvBuff[j];
    }
}
void Boundary::interiorYm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    int dest = mpi_grid.get_ym_id();
    // Remember tag should be for a Yp send
    // std::cout << mpi_grid.get_id() << ": Ym Waiting data from: "<< dest << std::endl;
    MPI_Recv(ymRecvBuff, lx, MPI_DOUBLE, dest, 1231, MPI_COMM_WORLD, &bstatus);
    // std::cout << mpi_grid.get_id() << ": Recieved data from: "<< dest << std::endl;
    // Set A's ghost nodes
    for(int i=0; i<lx; i++){
        A[0][i+1] = ymRecvBuff[i];
    }
}
void Boundary::interiorXm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);

    int dest = mpi_grid.get_xm_id();
    // Remember tag should be for a Xp send
    // std::cout << mpi_grid.get_id() << ": Xm Waiting data from: "<< dest << std::endl;
    MPI_Recv(xmRecvBuff, ly, MPI_DOUBLE, dest, 1232, MPI_COMM_WORLD, &bstatus);
    // std::cout << mpi_grid.get_id() << ": Recieved data from: "<< dest << std::endl;
    // Set A's ghost nodes
    for(int j=0; j<ly; j++){
        A[j+1][0] = xmRecvBuff[j];
    }
}

/**
 * Dirichlet boundary conditions. Use bndValues array to
 * change the value at the boundary. Only supports constant
 * value for entire side.
 */
 void Boundary::dirichletYp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int i=0; i<=lx+1; i++){
        A[ly+1][i] = bndValues[0];
    }
}
void Boundary::dirichletXp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int j=0; j<=ly+1; j++){
        A[j][lx+1] = bndValues[1];
    }
}
void Boundary::dirichletYm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int i=0; i<=lx+1; i++){
        A[0][i] = bndValues[2];
    }
}
void Boundary::dirichletXm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int j=0; j<=ly+1; j++){
        A[j][0] = bndValues[3];
    }
}

/**
 * Open boundary conditions. These require set_discretization and
 * set_material_field to be called prior since we will need to approx.
 * a gradient.
 */
void Boundary::openYp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    // Remeber K does not have ghost nodes
    for(int i=0; i<lx; i++){
        A[ly+1][i+1] = A[ly+1][i+1] - K[ly-1][i]*dt*(A[ly+1][i+1] - A[ly][i+1])/dy;
    }
}
void Boundary::openXp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    // Remeber K does not have ghost nodes
    for(int j=0; j<ly; j++){
        A[j+1][lx+1] = A[j+1][lx+1] - K[j][lx-1]*dt*(A[j+1][lx+1] - A[j+1][lx])/dx;
    }
}
void Boundary::openYm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    // Remeber K does not have ghost nodes
    for(int i=0; i<lx; i++){
        A[0][i+1] = A[0][i+1] + K[0][i]*dt*(A[1][i+1] - A[0][i+1])/dy;
    }
}
void Boundary::openXm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    // Remeber K does not have ghost nodes
    for(int j=0; j<ly; j++){
        A[j+1][0] = A[j+1][0] + K[j][0]*dt*(A[j+1][1] - A[j+1][0])/dx;
    }
}

/**
 * Locally periodic boundary conditions. Used
 * when the neighboring domain is the current.
 * Good when using small processes, not used otherwise
 */
void Boundary::localPeriodicYp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int i=0; i<=lx+1; i++){
        A[ly+1][i] = A[1][i];
    }
}
void Boundary::localPeriodicXp(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int j=0; j<=ly+1; j++){
        A[j][lx+1] = A[j][1];
    }
}
void Boundary::localPeriodicYm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int i=0; i<=lx+1; i++){
        A[0][i] = A[ly+1][i];
    }
}
void Boundary::localPeriodicXm(double**& A){
    int lx, ly;
    mpi_grid.get_domain_size(&lx, &ly);
    // Set ghost nodes
    for(int j=0; j<=ly+1; j++){
        A[j][0] = A[j][lx+1];
    }
}