#include "data_io.h"

DataIO::DataIO(MPI_Grid mpi_grid0){
    mpi_grid = mpi_grid0;
    mpi_grid.get_global_domain(&nx, &ny);

    nProc = mpi_grid.get_nproc();
}

/**
 * Collects state variable data from all local tasks and
 * writes to binary file
 * @param A - 2D array to be saved
 * @param t0 - double variable of current time-step value
 */
void DataIO::save_state(double**& A, double t0){

    if(mpi_grid.get_id() == 0){
        double** MA;
        alloc_2d_array(ny, nx, MA);

        int lSize[2];
        int lOrigin[2];
        mpi_grid.get_proc_domain(0, lOrigin, lSize);
        // First add local task to global array
        for(int i=0; i<lSize[0]; i++){
            for(int j=0; j<lSize[1]; j++){
                MA[i+lOrigin[0]][j+lOrigin[1]] = A[i+1][j+1];
            }
        }

        // Now recieve information from all other processes.
        for(int p=1; p<nProc; p++){
            // This returns matrix size so arrays are in [y,x] order
            mpi_grid.get_proc_domain(p, lOrigin, lSize);
            
            double** recvBuff;
            alloc_2d_array(lSize[0], lSize[1], recvBuff);
            MPI_Recv(&recvBuff[0][0], lSize[0]*lSize[1], MPI_DOUBLE, p, 1235, MPI_COMM_WORLD, &dstatus);
            // Recv buff to global array
            for(int i=0; i<lSize[0]; i++){
                for(int j=0; j<lSize[1]; j++){
                    MA[i+lOrigin[0]][j+lOrigin[1]] = recvBuff[i][j];
                }
            }
            // Delete recv buffer
            delete [] recvBuff[0];
            delete [] recvBuff;
        }

        // With all data collected, we can save the global to binary
        // First set file name with time-step
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << t0;
        std::string mystring = "data/u"+ss.str()+".dat";
        const char *file_name = mystring.c_str();
        save_array_binary(MA, nx, ny, file_name);

        mystring = "data/u"+ss.str()+".csv";
        const char *file_name2 = mystring.c_str();
        // save_array_text(MA, nx, ny, file_name2);

        // Delete global array
        delete [] MA[0];
        delete [] MA;

    }else{
        // If not master task send local information
        int lSize[2];
        int lOrigin[2];
        int myId = mpi_grid.get_id();
        mpi_grid.get_proc_domain(myId, lOrigin, lSize);

        double** sendBuff;
        alloc_2d_array(lSize[0], lSize[1], sendBuff);
        for(int i=0; i<lSize[0]; i++){
            for(int j=0; j<lSize[1]; j++){
                // Do not send ghost nodes
                sendBuff[i][j] = A[i+1][j+1];
            }
        }
        // Send data off to master
        MPI_Send(&sendBuff[0][0], lSize[0]*lSize[1], MPI_DOUBLE, 0, 1235, MPI_COMM_WORLD);

        // Delete send buffer
        delete [] sendBuff[0];
        delete [] sendBuff;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * Loads hetergenous material field from binary file
 * @param K - 2D material array to be loaded
 * @param fileName - character array of file name to load
 */
void DataIO::load_material(double**& K, const char* fileName){
    // Have only proc 0 read from file and then distribute
    if(mpi_grid.get_id() == 0){
        double** MK;
        alloc_2d_array(ny, nx, MK);

        // Read file from disk
        load_array_binary(MK, nx, ny, fileName);

        // Transfer to current tasks material array first
        int lSize[2];
        int lOrigin[2];
        mpi_grid.get_proc_domain(0, lOrigin, lSize);

        // First add local task to global array
        for(int i=0; i<lSize[0]; i++){
            for(int j=0; j<lSize[1]; j++){
                // No ghost nodes in K
                K[i][j] = MK[i][j];
            }
        }

        // Now send information from all other processes.
        // Use blocking communication here, not ideal but
        // only called once so not worth over complicating.
        for(int p=1; p<nProc; p++){
            // This returns matrix size so arrays are in [y,x] order
            mpi_grid.get_proc_domain(p, lOrigin, lSize);
            
            double** sendBuff;
            alloc_2d_array(lSize[0], lSize[1], sendBuff);
            for(int i=0; i<lSize[0]; i++){
                for(int j=0; j<lSize[1]; j++){
                    // Get data from global field
                    sendBuff[i][j] = MK[i+lOrigin[0]][j+lOrigin[1]];
                }
            }
            MPI_Send(&sendBuff[0][0], lSize[0]*lSize[1], MPI_DOUBLE, p, 1235, MPI_COMM_WORLD);

            // Delete recv buffer
            delete [] sendBuff[0];
            delete [] sendBuff;
        }
        // Delete global array
        delete [] MK[0];
        delete [] MK;
    }else{
        // If not master task recieve local material field
        int lSize[2];
        int lOrigin[2];
        int myId = mpi_grid.get_id();
        mpi_grid.get_proc_domain(myId, lOrigin, lSize);

        double** recvBuff;
        alloc_2d_array(lSize[0], lSize[1], recvBuff);
        
        // Recv data from master
        MPI_Recv(&recvBuff[0][0], lSize[0]*lSize[1], MPI_DOUBLE, 0, 1235, MPI_COMM_WORLD, &dstatus);
        for(int i=0; i<lSize[0]; i++){
            for(int j=0; j<lSize[1]; j++){
                // No ghost nodes in K
                K[i][j] = recvBuff[i][j];
            }
        }
        // Delete send buffer
        delete [] recvBuff[0];
        delete [] recvBuff;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * Saves 2D array to binary file.
 * https://stackoverflow.com/questions/16002680/writing-reading-2d-array-to-binary-file-c
 * @param A - 2D array to be saved
 * @param nx0, ny0 - integers of column and row dims respectively
 * @param fileName - character array of file name
 */
void DataIO::save_array_binary(double**& A, int nx0, int ny0, const char* fileName){

    std::cout << "Saving binary file: "<< fileName << std::endl;

    boost::filesystem::path dir(fileName);
    if(boost::filesystem::create_directories(dir.parent_path())) {
        std::cout << "Created folder: "<< dir.parent_path().c_str() << std::endl;
    }

    ofp.open(fileName, std::ios::out | std::ios::binary);
    if (!ofp){
        std::cout << "Error opening binary file." << std::endl;
        std::cout << "Issues with file: " << fileName << std::endl;
    }else{
        ofp.write(reinterpret_cast<char *>(&A[0][0]), (nx0*ny0)*sizeof(double));
        ofp.close();
    }
}

/**
 * Loads 2D array from binary file produced from matlab.
 * https://stackoverflow.com/questions/16002680/writing-reading-2d-array-to-binary-file-c
 * @param A - 2D array to be saved
 * @param nx0, ny0 - integers of column and row dims respectively
 * @param fileName - character array of file name
 */
void DataIO::load_array_binary(double**& A, int nx0, int ny0, const char* fileName){

    std::cout << "Reading binary file: "<< fileName << std::endl;
    FILE* fid;
    fid = fopen(fileName,"rb");
    if (!ifp){
        std::cout << "Error opening binary file." << std::endl;
        std::cout << "Issues with file: " << fileName << std::endl;
    }else{
        std::cout << nx0 << ny0 <<std::endl;
        fread(reinterpret_cast<char *>(&A[0][0]), (nx0*ny0)*sizeof(double), 1, fid);
        ifp.close();
    }
}

/**
 * Saves 2D array to text file. Only recommended for debugging.
 * During simulation, use binary format.
 * @param A - 2D array to be saved
 * @param nx0, ny0 - integers of column and row dims respectively
 * @param fileName - character array of file name
 */
void DataIO::save_array_text(double**& A, int nx0, int ny0, const char* fileName){

    std::cout << "Saving text file: "<< fileName << std::endl;
    ofp.open(fileName, std::ios::out);
    if (!ofp){
        std::cout << "Error opening text file." << std::endl;
        std::cout << "File: " << fileName << std::endl;
    }else{
        // Write array to file
        for(int i=0; i<ny0; i++){
            for(int j=0; j<nx0; j++){
                ofp << A[i][j] << " ";
            }
            ofp << std::endl;
        }
        ofp.close();
    }
}