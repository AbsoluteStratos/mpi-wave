#include  "common.h"

int alloc_2d_array(
       int   row,
       int   col,
       double**& pptr)
{
       int    i;

       pptr = new (std::nothrow) double*[row]; 
       if(nullptr == pptr){
           std::cout<<"Error in alloc_2d_array()\n"; 
	   return ERROR; 
       }

       pptr[0] = new (std::nothrow) double [row*col]; 
       if(nullptr == pptr[0]){
           std::cout<<"Error in alloc_2d_array() 2\n"; 
	   delete [] pptr; 
	   return ERROR; 
       }

       for(i = 1; i < row; i++)
       {
           pptr[i] = pptr[i-1]+col; 
       }

       return SUCCEED; 
}

int alloc_2d_int_array(
       int   row,
       int   col,
       int**& pptr)
{
       int    i;

       pptr = new (std::nothrow) int*[row]; 
       if(nullptr == pptr){
           std::cout<<"Error in alloc_2d_array()\n"; 
	   return ERROR; 
       }

       pptr[0] = new (std::nothrow) int [row*col]; 
       if(nullptr == pptr[0]){
           std::cout<<"Error in alloc_2d_array() 2\n"; 
	   delete [] pptr; 
	   return ERROR; 
       }

       for(i = 1; i < row; i++)
       {
           pptr[i] = pptr[i-1]+col; 
       }

       return SUCCEED; 
}

int alloc_1d_array(
	int    size,
	double* &ptr)
{
       ptr = new (std::nothrow) double [size]; 
       if(nullptr == ptr){
           std::cout<<"Error in alloc_1d_array()\n"; 
	   return ERROR; 
       }
       return SUCCEED; 
}


void create_mixed_xfer_arrays(
	int     id,
        int     p,
        int     n,
        int     *&count,
        int     *&disp)
{
	int     i;

        // *count = (int*) my_malloc(id, p*sizeof(int));
        // *disp = (int*) my_malloc(id, p*sizeof(int));

        count = new int[p]; 
        disp = new int[p]; 

        (count)[0] = BLOCK_SIZE(0,p,n);
        (disp)[0] = 0;

        for(i = 1; i < p; i++)
        {
            (disp)[i] = (disp)[i-1] + (count)[i-1];
            (count)[i] = BLOCK_SIZE(i,p,n);
        }
}


