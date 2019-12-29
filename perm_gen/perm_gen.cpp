#include<iostream>
#include<fstream>
#include <sstream>
#include<stdio.h> 
#include <string> 
#include <cstring>
using namespace std;
int main()
{
   union 
   { 
	   char b[8]; 
	   double d; 
   };
   //open the file
   ifstream rf("./mat.bin", ios::out | ios::binary);
   if(!rf) 
   {
      cout << "Cannot open file!" << endl;
      return 1;
   }
   
   // read the full size of the file (nx*ny*size of each value (for int =2, float =4, double =8))
   char A1[32768]; //char here
   rf.read(A1, (4096)*sizeof(double));//reinterpret_cast<char *>(&A[0]), (4096)*sizeof(double));
   int A[32768]; //to convert to int in bytes
   for (int i=0;i<32768;i++)
   {
		   A[i] = (int)A1[i];
		   A[i] = 255 + A[i] +1; //convert to int 255+a(i)+1
    }
   
   double A2[4096]; //in double values
   for (int j=0;j<4096;j++)
   {
	   int e[8]; //each value has 8 bytes
	   for (int k=0;k < 8;k++) //each byte loop
	   {
	   	e[k] = A[(j*8)+k]; //value in int (max =255)
		   b[k] = 0xff & (unsigned int) e[k]; //write to hex format
	   }
	   A2[j] = d; //in double value using union definition
   }
   //print function of A2 decimal double value here
   for (int i =0;i<4096;i++)
   {
	   cout << A2[i] << " "; //just for display
   }
   //issue with reading -> report here
   if(!rf.good())
   {
      cout << "Error occurred at reading time!" << endl;
      return 1;
   }
   //close file
   rf.close();
   return 0;
}
