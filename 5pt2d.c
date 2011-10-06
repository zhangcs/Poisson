/*
 *  This is a program to generate linear system for 2D Possion Problem.
 *
 *     Consider a two-dimensional Possion equation
 *
 *          / -u_{xx}-u_{yy} = f(x,y)     in \Omega = (0,1)X(0,1)
 *          \  u = 0                      on \partial\Omega
 *
 *  where f(x,y = 2*pi^2*sin(pi*x)*sin(pi*y), 
 *  and the solution function can be expressed by
 * 
 *             u(x,y) = sin(pi*x)*sin(pi*y)
 *
 *  Created by peghoty 2010/08/04
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 */

#include "fsls.h"

int 
main( int argc, char *argv[])
{
   struct timeval tStart,tEnd;
   int TTest       = 1;
   int arg_index   = 0;
   int print_usage = 0;  

   char *MatFile = NULL;
   char *RhsFile = NULL;
   char *SolFile = NULL;
   char filename[120];  
   
   fsls_BandMatrix *A    = NULL;
   fsls_CSRMatrix  *Acsr = NULL;
   fsls_XVector    *b    = NULL;
   fsls_XVector    *u    = NULL;
   
   int nx,ny,ngrid;
   
   nx = 1024;
   ny = 1024; 
 
   while (arg_index < argc)
   {
       if ( strcmp(argv[arg_index], "-nx") == 0 )
       {
           arg_index ++;
           nx = atoi(argv[arg_index++]);
       }
       if ( strcmp(argv[arg_index], "-ny") == 0 )
       {
           arg_index ++;
           ny = atoi(argv[arg_index++]);
       }                                                                              
       else if ( strcmp(argv[arg_index], "-help") == 0 )
       {
           print_usage = 1;
           break;
       }
       else
       {
           arg_index ++;
       }         
   } 
 
   if (print_usage)
   {
         printf("\n");
         printf("  Usage: %s [<options>]\n", argv[0]);
         printf("\n");       
         printf("  -nx     <val> : number of interier nodes in x-direction[default: 10]\n");
         printf("  -ny     <val> : number of interier nodes in y-direction[default: 10]\n");       
         printf("  -help         : using help message\n\n");
         exit(1);
   }
   
   ngrid = nx*ny;
   printf("\n +++++++++++++ (nx,ny) = (%d,%d)  ngrid = %d +++++++++++\n\n",nx,ny,ngrid);

   MatFile = "../input/mat_";
   RhsFile = "../input/rhs_";
   SolFile = "../input/sol_"; 

  /*-----------------------------------------------------
   * construct a linear system
   *----------------------------------------------------*/  
   if (TTest) GetTime(tStart); 
     
   fsls_BuildLinearSystem_5pt2d(nx, ny, &A, &b, &u);
   
   if (TTest) 
   {
      GetTime(tEnd);
      printf("\n >>> total time: %.3f seconds\n\n",mytime(tStart,tEnd));       
   }   

   sprintf(filename, "%s%dX%d",MatFile,nx,ny);
   fsls_Band2CSRMatrix(A, &Acsr);
   fsls_CSRMatrixPrint(Acsr,filename);
   
   sprintf(filename, "%s%dX%d",RhsFile,nx,ny);
   fsls_XVectorPrint(b, filename);
   
   sprintf(filename, "%s%dX%d",SolFile,nx,ny);
   fsls_XVectorPrint(u, filename);
   
  /*------------------------------------------------------
   * added for Chensong's SAMG-testing
   *-----------------------------------------------------*/
   fsls_WriteSAMGData(Acsr, b, u);
   
  /*------------------------------------------------------
   * free some staff
   *-----------------------------------------------------*/   
   fsls_BandMatrixDestroy(A);
   fsls_CSRMatrixDestroy(Acsr);
   fsls_XVectorDestroy(b);
   fsls_XVectorDestroy(u);
   
   return(0);
}



