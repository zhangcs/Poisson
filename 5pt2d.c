/*
 *  This is a program to generate linear system for 2D Possion Problem.
 *
 *     Consider a two-dimensional Possion equation
 *
 *          /du/dt-u_{xx}-u_{yy} = f(x,y,t) in [0,1] X \Omega = (0,1)X(0,1)
 *          |          u_(x,y,0) = 0        in \Omega
 *          \                  u = 0        on (0,1] X \partial\Omega
 *
 *  where f(x,y,t) = 2*pi^2*sin(pi*x)*sin(pi*y)*t + sin(pi*x)*sin(pi*y),
 *  and the solution function can be expressed by
 * 
 *             u(x,y,t) = sin(pi*x)*sin(pi*y)*t
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
   
   int nx,ny,ngrid,nt,i,j;
   
   nx = 10;
   ny = 10;
	 nt = 20;
	 double dt = 1./nt;
 
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
       if ( strcmp(argv[arg_index], "-nt") == 0 )
       {
           arg_index ++;
           nt = atoi(argv[arg_index++]);
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
         printf("  -nt     <val> : number of interier nodes in t-direction[default: 20]\n");
         printf("  -help         : using help message\n\n");
         exit(1);
   }
   
   ngrid = nx*ny;
   printf("\n +++++++++++++ (nx,ny,nt) = (%d,%d,%d)  ngrid = %d +++++++++++\n\n",nx,ny,nt,ngrid);

   MatFile = "./mat_";
   RhsFile = "./rhs_";
   SolFile = "./sol_"; 

  /*-----------------------------------------------------
   * construct a linear system
   *----------------------------------------------------*/  
   if (TTest) GetTime(tStart); 
     
   fsls_BuildLinearSystem_5pt2d(nt, nx, ny, &A, &b, &u);
   
   if (TTest) 
   {
      GetTime(tEnd);
      printf("\n >>> total time: %.3f seconds\n\n",mytime(tStart,tEnd));       
   }   

   sprintf(filename, "%s%dX%d.dat",MatFile,nx,ny);
   fsls_Band2CSRMatrix(A, &Acsr);
	 //newly added 2011/12/11 by feiteng,use lapack routine as solver
	 double *A_full = NULL, *B = NULL;
	 int LDA, LDB, NRHS = 1;
	 int *IPIV, INFO[1];
	 fsls_CSR2FullMatrix(Acsr, &A_full);
	 fsls_dtMatrix(dt, ngrid, ngrid, A_full);

	 IPIV = ( int *) malloc ( sizeof (int) * ngrid ) ;
	 if( IPIV == NULL ) {
		 printf("Allocation Error\n");
		 exit(1);
    }
	 LDA = LDB = ngrid;
	 B = ( double *) malloc ( sizeof(double) * ngrid);
	 if( B == NULL ){
		 printf("Allocation Error\n");
		 exit(1);
	 }
//	 memcpy(B, fsls_XVectorData(b), sizeof(double) * ngrid);
	 memset(B, 0X0, ngrid*sizeof(double));
	 dgetrf_(&ngrid, &ngrid, A_full, &LDA, IPIV, INFO);
	 for (i = 0;i < nt; ++i)
	 {
		 double err = 0.;
		 double u_ = 0.;
		 for (j = 0;j < ngrid; ++j)
		 {
			 B[j] += b->data[i*ngrid+j]*dt;
		 }
		 dgetrs_("N", &ngrid, &NRHS, A_full, &LDA, IPIV, B, &LDB, INFO);
		 for (j = 0;j < ngrid; ++j)
		 {
			 err += (B[j]-u->data[i*ngrid+j])*(B[j]-u->data[i*ngrid+j]);
			 u_ += (u->data[i*ngrid+j])*(u->data[i*ngrid+j]);
		 }
		 printf("...%f\n",err/u_);
	 }

   fsls_CSRMatrixPrint(Acsr,filename);
   
   sprintf(filename, "%s%dX%d.dat",RhsFile,nx,ny);
   fsls_XVectorPrint(b, filename);
   
   sprintf(filename, "%s%dX%d.dat",SolFile,nx,ny);
   fsls_XVectorPrint(u, filename);
   
  /*------------------------------------------------------
   * added for Chensong's SAMG-testing
   *-----------------------------------------------------*/
   // fsls_WriteSAMGData(Acsr, b, u);
   
  /*------------------------------------------------------
   * free some staff
   *-----------------------------------------------------*/   
   fsls_BandMatrixDestroy(A);
   fsls_CSRMatrixDestroy(Acsr);
   fsls_XVectorDestroy(b);
   fsls_XVectorDestroy(u);

	 free(A_full);
	 free(B);
	 free(IPIV);
  
   return(0);
}



