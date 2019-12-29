#if !defined(_COMM_H)
#define _COMM_H

#include <iostream>
#include <fstream> 
#include <new>
#include <cmath>
#include <cstdlib>
#include <cstring>

#define PI        3.14159265
#define ZERO      1.0E-20
#define SUCCEED   1
#define ERROR     0

// #define N      4  /*Number of equations*/

#if !defined(__cplusplus)
#define true   1
#define false  0
#endif

int alloc_2d_array(int,int,double**&); 
int alloc_2d_int_array(int,int,int**&); 
int alloc_1d_array(int,double*&); 
void create_mixed_xfer_arrays(int,int,int,int*&,int*&); 

#define DATA_MSG     0
#define PROMPT_MSG   1
#define RESPONSE_MSG   2

#define OPEN_FILE_ERROR   -1
#define MALLOC_ERROR      -2
#define TYPE_ERROR        -3

#define  BLOCK_LOW(id, p, n)     ((id)*(n)/(p))
#define  BLOCK_HIGH(id, p, n)    (BLOCK_LOW((id)+1,p,n)-1)
#define  BLOCK_SIZE(id, p, n)    (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n)+1)
#define  BLOCK_OWNER(j, p, n)    (((p)*((j)+1)-1)/(n))
#define  PTR_SIZE                (sizeof(void*))
#define  CEILING(i,j)            (((i)+(j)-1)/(j))

#endif 